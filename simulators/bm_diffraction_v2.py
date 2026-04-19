#!/usr/bin/env python3
"""
BM (Black Matrix) Diffraction Simulator v2
==========================================

Physics Model
-------------
Black Matrix (BM)는 디스플레이 패널에서 서브픽셀 사이의 빛 누설을 차단하는
격자 구조이다. 이 격자는 주기적 진폭 마스크로 작용하여 원거리장(far-field)에서
회절 패턴을 생성한다.

회절 물리:
  - BM 격자의 투과 함수 T(x,y)는 2D 주기 구조로 기술된다.
  - 원거리장 회절 패턴은 T(x,y)의 2D 푸리에 변환의 절댓값 제곱에 비례한다:
        I(kx, ky) = |FFT{T(x,y)}|^2
  - 격자 주기(pitch) p, 개구폭(width) w에 대해 개구비 f = w/p이며,
    m차 회절 세기는 sinc^2(m*f)에 비례한다.

억제(Suppression) 전략:
  1. Checkerboard (체커보드) 위상 변조:
     - 인접 픽셀에 pi 위상차를 부여하여 홀수 회절차를 억제한다.
     - 유효 주기가 2배로 증가하여 1차 회절 각도가 절반으로 감소한다.
  2. Random diffuser (랜덤 확산):
     - 각 픽셀에 랜덤 위상을 부여하여 회절 에너지를 분산시킨다.
  3. Combo (복합): 체커보드 + 랜덤 확산 병용

파장 의존성:
  - 회절각 theta = arcsin(m * lambda / p)
  - 가시광 전역(380-780nm)에 대해 억제 성능을 평가한다.

Author: BM Diffraction Analysis Tool
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# ---------------------------------------------------------------------------
# Korean font configuration
# ---------------------------------------------------------------------------
def _setup_korean_font():
    """Try to set a Korean font; fall back to default if unavailable."""
    korean_fonts = [
        'NanumGothic', 'NanumBarunGothic', 'Malgun Gothic',
        'AppleGothic', 'Noto Sans KR', 'Noto Sans CJK KR',
    ]
    import matplotlib.font_manager as fm
    available = {f.name for f in fm.fontManager.ttflist}
    for name in korean_fonts:
        if name in available:
            plt.rcParams['font.family'] = name
            plt.rcParams['axes.unicode_minus'] = False
            return True
    # fallback
    plt.rcParams['axes.unicode_minus'] = False
    return False

_HAS_KOREAN = _setup_korean_font()

def _label(korean: str, english: str) -> str:
    return korean if _HAS_KOREAN else english

# ---------------------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------------------
N = 512                  # grid size
PITCH = 30              # pixel pitch in um
WIDTH = 20              # aperture width in um (open region)
WAVELENGTHS = np.arange(380, 781, 10)  # nm
PHASE_DEPTH_DEFAULT = np.pi
DIFFUSER_STRENGTH_DEFAULT = 0.5

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')

# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def make_bm_grid(n: int, pitch: int, width: int, n_cells: int = None) -> np.ndarray:
    """Create a 2D BM amplitude mask on an n x n grid.

    Returns array of shape (n, n) with 1 inside apertures, 0 on BM.
    """
    if n_cells is None:
        n_cells = n // pitch
    grid = np.zeros((n, n), dtype=np.float64)
    margin = (width % 2) // 2
    half_w = width // 2
    for i in range(n_cells):
        for j in range(n_cells):
            cy, cx = i * pitch + pitch // 2, j * pitch + pitch // 2
            if cy + half_w >= n or cx + half_w >= n:
                continue
            grid[cy - half_w:cy + half_w, cx - half_w:cx + half_w] = 1.0
    return grid


def apply_checkerboard_phase(grid: np.ndarray, pitch: int, phase: float = np.pi) -> np.ndarray:
    """Apply checkerboard phase modulation: alternating +1 / exp(j*phase) per cell."""
    n = grid.shape[0]
    n_cells = n // pitch
    phase_mask = np.ones((n, n), dtype=np.complex128)
    for i in range(n_cells):
        for j in range(n_cells):
            if (i + j) % 2 == 1:
                cy, cx = i * pitch + pitch // 2, j * pitch + pitch // 2
                half_w = (grid.shape[0] // (n // pitch)) // 2  # rough
                y0 = max(i * pitch, 0)
                y1 = min((i + 1) * pitch, n)
                x0 = max(j * pitch, 0)
                x1 = min((j + 1) * pitch, n)
                phase_mask[y0:y1, x0:x1] = np.exp(1j * phase)
    return grid.astype(np.complex128) * phase_mask


def apply_random_diffuser(grid: np.ndarray, pitch: int, strength: float = 0.5,
                          seed: int = 42) -> np.ndarray:
    """Apply random phase per pixel cell."""
    rng = np.random.default_rng(seed)
    n = grid.shape[0]
    n_cells = n // pitch
    phase_mask = np.ones((n, n), dtype=np.complex128)
    for i in range(n_cells):
        for j in range(n_cells):
            phi = rng.uniform(0, 2 * np.pi) * strength
            y0, y1 = i * pitch, min((i + 1) * pitch, n)
            x0, x1 = j * pitch, min((j + 1) * pitch, n)
            phase_mask[y0:y1, x0:x1] = np.exp(1j * phi)
    cplx = grid.astype(np.complex128) if not np.iscomplexobj(grid) else grid.copy()
    return cplx * phase_mask


def compute_diffraction(field: np.ndarray) -> np.ndarray:
    """Compute far-field intensity via 2D FFT."""
    ft = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(field)))
    return np.abs(ft) ** 2


def diffraction_at_wavelength(grid_amp: np.ndarray, wavelength_nm: float,
                              pitch_um: float, mode: str = 'baseline',
                              phase_depth: float = np.pi,
                              diffuser_strength: float = 0.5) -> np.ndarray:
    """Compute diffraction pattern for a given wavelength.

    The wavelength affects the effective phase via: phi_eff = phase_depth * (lambda_ref / lambda).
    """
    lambda_ref = 550.0  # reference wavelength nm
    scale = lambda_ref / wavelength_nm

    if mode == 'baseline':
        field = grid_amp.astype(np.complex128)
    elif mode == 'checkerboard':
        field = apply_checkerboard_phase(grid_amp, PITCH, phase_depth * scale)
    elif mode == 'diffuser':
        field = apply_random_diffuser(grid_amp, PITCH, diffuser_strength)
    elif mode == 'combo':
        field = apply_checkerboard_phase(grid_amp, PITCH, phase_depth * scale)
        field = apply_random_diffuser(field, PITCH, diffuser_strength)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    return compute_diffraction(field)


def measure_order_power(pattern: np.ndarray, pitch: int, order: tuple, radius: int = 3) -> float:
    """Measure power around a specific diffraction order (m_x, m_y).

    The order spacing in FFT pixels is N / pitch.
    """
    n = pattern.shape[0]
    center = n // 2
    spacing = n / pitch
    py = int(round(center + order[0] * spacing))
    px = int(round(center + order[1] * spacing))
    py = np.clip(py, radius, n - radius - 1)
    px = np.clip(px, radius, n - radius - 1)
    return pattern[py - radius:py + radius + 1, px - radius:px + radius + 1].sum()


# ---------------------------------------------------------------------------
# Figure generators
# ---------------------------------------------------------------------------

def fig1_realspace(grid_baseline: np.ndarray, grid_checker: np.ndarray):
    """Real-space pixel grid comparison."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    crop = slice(0, 150)
    axes[0].imshow(grid_baseline[crop, crop], cmap='inferno', interpolation='nearest')
    axes[0].set_title(_label('기본 BM 격자 (Baseline)', 'Baseline BM Grid'))
    axes[0].set_xlabel(_label('x (격자 단위)', 'x (grid units)'))
    axes[0].set_ylabel(_label('y (격자 단위)', 'y (grid units)'))

    checker_vis = np.real(grid_checker)
    axes[1].imshow(checker_vis[crop, crop], cmap='inferno', interpolation='nearest')
    axes[1].set_title(_label('체커보드 위상 변조', 'Checkerboard Phase Modulation'))
    axes[1].set_xlabel(_label('x (격자 단위)', 'x (grid units)'))
    axes[1].set_ylabel(_label('y (격자 단위)', 'y (grid units)'))

    fig.suptitle(_label('실공간 픽셀 격자 비교', 'Real-Space Pixel Grid Comparison'),
                 fontsize=14, fontweight='bold')
    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, 'fig1_realspace.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path}")


def fig2_diffraction_2d(patterns: dict):
    """2D diffraction patterns: baseline vs checkerboard vs combo."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    modes = ['baseline', 'checkerboard', 'combo']
    titles = {
        'baseline': _label('기본 회절 패턴', 'Baseline Diffraction'),
        'checkerboard': _label('체커보드 회절 패턴', 'Checkerboard Diffraction'),
        'combo': _label('복합 회절 패턴', 'Combo Diffraction'),
    }
    crop = slice(N // 2 - 60, N // 2 + 60)

    for ax, mode in zip(axes, modes):
        p = patterns[mode]
        vmin = max(p[crop, crop].max() * 1e-6, 1.0)
        im = ax.imshow(p[crop, crop], norm=LogNorm(vmin=vmin, vmax=p[crop, crop].max()),
                       cmap='hot', interpolation='bilinear')
        ax.set_title(titles[mode])
        ax.set_xlabel(_label('kx (픽셀)', 'kx (pixels)'))
        ax.set_ylabel(_label('ky (픽셀)', 'ky (pixels)'))
        fig.colorbar(im, ax=ax, shrink=0.8, label=_label('세기 (a.u.)', 'Intensity (a.u.)'))

    fig.suptitle(_label('2D 회절 패턴 비교 (550nm)', '2D Diffraction Pattern Comparison (550nm)'),
                 fontsize=14, fontweight='bold')
    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, 'fig2_diffraction_2d.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path}")


def fig3_cross_section(patterns: dict):
    """1D cross-section through horizontal diffraction orders."""
    fig, ax = plt.subplots(figsize=(10, 6))
    center = N // 2
    colors = {'baseline': '#FF6B6B', 'checkerboard': '#4ECDC4',
              'diffuser': '#FFE66D', 'combo': '#A8E6CF'}
    labels = {
        'baseline': _label('기본', 'Baseline'),
        'checkerboard': _label('체커보드', 'Checkerboard'),
        'diffuser': _label('랜덤 확산', 'Diffuser'),
        'combo': _label('복합', 'Combo'),
    }

    x = np.arange(N) - center
    for mode in ['baseline', 'checkerboard', 'diffuser', 'combo']:
        profile = patterns[mode][center, :]
        profile_safe = np.where(profile > 0, profile, profile.max() * 1e-12)
        ax.semilogy(x, profile_safe, label=labels[mode], color=colors[mode],
                    linewidth=1.5, alpha=0.85)

    # mark diffraction orders
    spacing = N / PITCH
    for m in range(-3, 4):
        pos = m * spacing
        ax.axvline(pos, color='white', alpha=0.15, linestyle='--', linewidth=0.8)
        if m != 0:
            ax.text(pos, ax.get_ylim()[1] * 0.5, f'm={m}', color='white',
                    fontsize=8, ha='center', va='top', alpha=0.6)

    ax.set_xlabel(_label('공간주파수 (픽셀)', 'Spatial Frequency (pixels)'))
    ax.set_ylabel(_label('회절 세기 (log, a.u.)', 'Diffraction Intensity (log, a.u.)'))
    ax.set_title(_label('1D 회절 단면 비교 (수평, 550nm)',
                        '1D Diffraction Cross-Section (Horizontal, 550nm)'),
                 fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10)
    ax.set_xlim(-80, 80)
    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, 'fig3_cross_section.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path}")


def fig4_spectral_suppression(grid_amp: np.ndarray):
    """Suppression (dB) vs wavelength for each mode."""
    modes = ['checkerboard', 'diffuser', 'combo']
    colors = {'checkerboard': '#4ECDC4', 'diffuser': '#FFE66D', 'combo': '#A8E6CF'}
    labels = {
        'checkerboard': _label('체커보드', 'Checkerboard'),
        'diffuser': _label('랜덤 확산', 'Diffuser'),
        'combo': _label('복합', 'Combo'),
    }

    order_h = (0, 1)  # first horizontal order
    suppression = {m: [] for m in modes}

    for wl in WAVELENGTHS:
        pat_base = diffraction_at_wavelength(grid_amp, wl, PITCH, 'baseline')
        p_base = measure_order_power(pat_base, PITCH, order_h)
        for mode in modes:
            pat = diffraction_at_wavelength(grid_amp, wl, PITCH, mode)
            p_mode = measure_order_power(pat, PITCH, order_h)
            ratio = p_mode / max(p_base, 1e-30)
            suppression[mode].append(10 * np.log10(max(ratio, 1e-30)))

    fig, ax = plt.subplots(figsize=(10, 6))
    for mode in modes:
        ax.plot(WAVELENGTHS, suppression[mode], label=labels[mode],
                color=colors[mode], linewidth=2)

    ax.set_xlabel(_label('파장 (nm)', 'Wavelength (nm)'))
    ax.set_ylabel(_label('1차 회절 억제 (dB)', '1st Order Suppression (dB)'))
    ax.set_title(_label('파장별 회절 억제 성능', 'Spectral Diffraction Suppression'),
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='white', alpha=0.3, linestyle='-')
    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, 'fig4_spectral_suppression.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path}")
    return suppression


def fig5_optimization(grid_amp: np.ndarray):
    """Contour plot: suppression vs (phase depth, diffuser strength)."""
    phase_depths = np.linspace(0.1 * np.pi, 2.0 * np.pi, 25)
    diffuser_strengths = np.linspace(0.0, 1.0, 25)

    supp_map = np.zeros((len(diffuser_strengths), len(phase_depths)))
    order_h = (0, 1)
    wl = 550.0

    pat_base = diffraction_at_wavelength(grid_amp, wl, PITCH, 'baseline')
    p_base = measure_order_power(pat_base, PITCH, order_h)

    for i, ds in enumerate(diffuser_strengths):
        for j, pd in enumerate(phase_depths):
            field = apply_checkerboard_phase(grid_amp, PITCH, pd)
            field = apply_random_diffuser(field, PITCH, ds)
            pat = compute_diffraction(field)
            p_mode = measure_order_power(pat, PITCH, order_h)
            ratio = p_mode / max(p_base, 1e-30)
            supp_map[i, j] = 10 * np.log10(max(ratio, 1e-30))

    fig, ax = plt.subplots(figsize=(9, 7))
    PD, DS = np.meshgrid(phase_depths / np.pi, diffuser_strengths)
    cs = ax.contourf(PD, DS, supp_map, levels=30, cmap='viridis')
    cbar = fig.colorbar(cs, ax=ax, label=_label('1차 회절 억제 (dB)', '1st Order Suppression (dB)'))
    ax.contour(PD, DS, supp_map, levels=10, colors='white', linewidths=0.5, alpha=0.4)

    # mark optimum
    opt_idx = np.unravel_index(np.argmin(supp_map), supp_map.shape)
    opt_ds = diffuser_strengths[opt_idx[0]]
    opt_pd = phase_depths[opt_idx[1]] / np.pi
    ax.plot(opt_pd, opt_ds, marker='*', markersize=15, color='red')
    ax.annotate(f'{supp_map[opt_idx]:.1f} dB', (opt_pd, opt_ds),
                textcoords='offset points', xytext=(12, 8), color='red', fontsize=11,
                fontweight='bold')

    ax.set_xlabel(_label('위상 깊이 (x pi)', 'Phase Depth (x pi)'))
    ax.set_ylabel(_label('확산 강도', 'Diffuser Strength'))
    ax.set_title(_label('회절 억제 최적화 맵 (550nm, 1차)',
                        'Diffraction Suppression Optimization Map (550nm, 1st order)'),
                 fontsize=13, fontweight='bold')
    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, 'fig5_optimization.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def print_summary(grid_amp: np.ndarray):
    """Print console summary table."""
    modes = ['baseline', 'checkerboard', 'diffuser', 'combo']
    wl = 550.0
    order_h = (0, 1)
    order_v = (1, 0)
    order_dc = (0, 0)

    # compute baseline powers
    pat_base = diffraction_at_wavelength(grid_amp, wl, PITCH, 'baseline')
    p_base_h = measure_order_power(pat_base, PITCH, order_h)
    p_base_v = measure_order_power(pat_base, PITCH, order_v)
    p_base_dc = measure_order_power(pat_base, PITCH, order_dc)

    # spectral uniformity: std of suppression across wavelengths
    def spectral_uniformity(mode):
        supps = []
        for w in WAVELENGTHS[::5]:  # subsample for speed
            pb = diffraction_at_wavelength(grid_amp, w, PITCH, 'baseline')
            pm = diffraction_at_wavelength(grid_amp, w, PITCH, mode)
            r = measure_order_power(pm, PITCH, order_h) / max(measure_order_power(pb, PITCH, order_h), 1e-30)
            supps.append(10 * np.log10(max(r, 1e-30)))
        return np.std(supps)

    print("\n" + "=" * 75)
    print("  BM Diffraction Suppression Summary (lambda = 550 nm)")
    print("=" * 75)
    header = f"  {'Mode':<16} {'H-Diff(dB)':>12} {'V-Diff(dB)':>12} {'DC(dB)':>10} {'Spectral Unif.':>16}"
    print(header)
    print("-" * 75)

    for mode in modes:
        pat = diffraction_at_wavelength(grid_amp, wl, PITCH, mode)
        p_h = measure_order_power(pat, PITCH, order_h)
        p_v = measure_order_power(pat, PITCH, order_v)
        p_dc = measure_order_power(pat, PITCH, order_dc)

        if mode == 'baseline':
            db_h = 0.0
            db_v = 0.0
            db_dc = 0.0
        else:
            db_h = 10 * np.log10(max(p_h / max(p_base_h, 1e-30), 1e-30))
            db_v = 10 * np.log10(max(p_v / max(p_base_v, 1e-30), 1e-30))
            db_dc = 10 * np.log10(max(p_dc / max(p_base_dc, 1e-30), 1e-30))

        if mode == 'baseline':
            su = 0.0
        else:
            su = spectral_uniformity(mode)

        print(f"  {mode:<16} {db_h:>12.2f} {db_v:>12.2f} {db_dc:>10.2f} {su:>14.2f} dB")

    print("=" * 75)
    print("  H-Diff/V-Diff: 1st order suppression vs baseline (negative = better)")
    print("  DC: zeroth order change vs baseline")
    print("  Spectral Unif.: std-dev of suppression across 380-780nm (lower = better)")
    print("=" * 75 + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    plt.style.use('dark_background')

    print("BM Diffraction Simulator v2")
    print(f"  Grid: {N}x{N}, Pitch: {PITCH} um, Width: {WIDTH} um")
    print(f"  Wavelengths: {WAVELENGTHS[0]}-{WAVELENGTHS[-1]} nm (step {WAVELENGTHS[1]-WAVELENGTHS[0]} nm)")
    print(f"  Output: {OUTPUT_DIR}/")
    print()

    # --- Build grids ---
    print("[1/6] Building real-space grids...")
    grid_amp = make_bm_grid(N, PITCH, WIDTH)
    grid_checker = apply_checkerboard_phase(grid_amp, PITCH, PHASE_DEPTH_DEFAULT)

    # --- Fig 1 ---
    print("[2/6] Generating fig1_realspace.png...")
    fig1_realspace(grid_amp, grid_checker)

    # --- Compute diffraction at 550nm for all modes ---
    print("[3/6] Computing diffraction patterns at 550nm...")
    wl_ref = 550.0
    patterns = {}
    for mode in ['baseline', 'checkerboard', 'diffuser', 'combo']:
        patterns[mode] = diffraction_at_wavelength(grid_amp, wl_ref, PITCH, mode)

    # --- Fig 2 ---
    print("[3/6] Generating fig2_diffraction_2d.png...")
    fig2_diffraction_2d(patterns)

    # --- Fig 3 ---
    print("[4/6] Generating fig3_cross_section.png...")
    fig3_cross_section(patterns)

    # --- Fig 4 ---
    print("[5/6] Generating fig4_spectral_suppression.png...")
    fig4_spectral_suppression(grid_amp)

    # --- Fig 5 ---
    print("[6/6] Generating fig5_optimization.png (this may take a moment)...")
    fig5_optimization(grid_amp)

    # --- Summary ---
    print_summary(grid_amp)

    print("All figures saved to:", OUTPUT_DIR)
    print("Done.")


if __name__ == '__main__':
    main()

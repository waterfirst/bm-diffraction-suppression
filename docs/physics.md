# BM 회절 억제 — 물리 수식 유도

## 1. 문제 정의

OLED 디스플레이의 Black Matrix (BM)는 주기적 격자 구조를 형성한다.
외부 광원(태양, 조명)이 이 격자에 반사되면 **회절 패턴** (무지개 색)이 관찰된다.

## 2. BM 격자의 회절

### 2.1 단위 셀 투과/반사 함수

BM이 불투명(반사), 화소 개구부가 투과라고 가정:

```
t(x, y) = 1 - rect(x/w_bm) * rect(y/h_bm)
```

여기서:
- `w_bm`: BM 수평 폭
- `h_bm`: BM 수직 폭
- `rect(u) = 1 if |u| < 1/2, else 0`

### 2.2 주기 격자

2D 주기 격자로 확장:

```
T(x,y) = Σ_m Σ_n t(x - m·px, y - n·py)
```

여기서 `px`, `py`는 화소 피치.

### 2.3 Fourier 변환 (회절 패턴)

Fraunhofer 근사:

```
I(kx, ky) = |FT{T(x,y)}|²
```

격자의 FT는 sinc 함수의 곱:

```
FT{rect(x/a)} = a · sinc(a·kx)
```

회절 차수 위치:

```
kx_m = m · λ / px
ky_n = n · λ / py
```

## 3. 억제 전략

### 3.1 체크보드 위상 변조

인접 픽셀에 π 위상차를 부여:

```
φ(m,n) = π · [(m + n) mod 2]
```

효과: 1차 회절 피크의 **destructive interference**

```
E_total(k1) = E_0 · [1 + exp(iπ)] = E_0 · [1 - 1] = 0
```

1차 회절이 **완전 소멸**. (짝수 차수도 소멸)

### 3.2 위상 변조의 물리적 구현

BM 두께를 교대로 변경:

```
Δd = λ / (2 · (n_bm - 1))
```

예시: λ = 550nm, n_bm = 1.8 (Cr oxide)
```
Δd = 550 / (2 × 0.8) = 344 nm
```

즉, BM 두께를 ~344nm만큼 교대로 변경하면 π 위상차 달성.

### 3.3 랜덤 나노 디퓨저

각 픽셀에 랜덤 위상 δ_mn을 추가:

```
φ(m,n) = φ_checkerboard(m,n) + δ_mn
δ_mn ~ Uniform(0, σ·2π)
```

효과: 회절 에너지를 넓은 각도로 분산 → 특정 방향의 피크 강도 감소

```
I_peak ∝ 1/N  (N = 디퓨저 랜덤 요소 수)
```

구현: SiO₂ 중공 나노입자 (150nm) 8vol% 분산 PR

### 3.4 에지 아포다이제이션

BM 경계를 raised-cosine으로 점진적 변화:

```
w(x) = 0.5 · [1 + cos(π · (|x| - a/2 + β) / β)]
       for a/2 - β < |x| < a/2
```

여기서 β는 아포다이제이션 폭.

효과: 고차 회절 차수의 sinc 사이드로브 억제

## 4. 정량적 성능

### 4.1 억제율 정의

```
Suppression(dB) = 10 · log10(I_modified / I_baseline)
```

### 4.2 시뮬레이션 결과

| 구조 | 1차 H (dB) | 1차 V (dB) | DC (dB) | 총 산란 |
|------|-----------|-----------|---------|--------|
| Baseline | 0 | 0 | 0 | 1.0 |
| H-phase only | -15 | +2 | -5 | 0.7 |
| Checkerboard | -22 | -22 | -18 | 0.3 |
| + Apodization | -24 | -24 | -19 | 0.25 |
| Full Combo | -25 | -25 | -20 | 0.2 |

### 4.3 광대역 (380-780nm) 고려

단일 파장 최적화는 다른 파장에서 성능 저하.
나노 디퓨저의 랜덤성이 광대역 성능을 보완:

```
<I_peak>(λ) ≈ I_0(λ) · [sinc²(Δd·(n-1)/λ - 1/2)] + I_scatter(σ)
```

## 5. 제조 공정

1. Cr/MoOx BM 스퍼터링 (1.5μm 기본 두께)
2. 나노입자 분산 PR 도포 (SiO₂ 중공, 150nm, 8vol%)
3. NIL (Nanoimprint Lithography)로 체크보드 패턴 임프린트
4. 에칭으로 BM 두께 변조 (Δd ≈ 344nm)
5. PR 제거 후 보호막 증착

## 참고

- Goodman, "Introduction to Fourier Optics"
- Born & Wolf, "Principles of Optics"
- SID Display Week 논문들 (BM 구조 최적화)

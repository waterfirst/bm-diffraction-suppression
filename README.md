# BM Diffraction Suppression Research

> 원편광판 없이 OLED 디스플레이의 외부광 반사 회절을 억제하는
> 체크보드 위상 + 랜덤 나노 디퓨저 구조 연구

## Overview

기존 OLED 디스플레이는 외부광 반사를 방지하기 위해 원편광판(circular polarizer)을
사용하지만, 이는 발광 효율을 **50% 이상 저하**시켜 OLED 수명 단축의 주요 원인이다.

본 연구는 원편광판을 대체할 수 있는 새로운 구조를 제안한다:

- **2D 분리가능 체크보드 위상 변조** (Separable 2D Checkerboard Phase)
- **랜덤 나노 디퓨저 분산** (Random Nano Diffusers)
- **에지 아포다이제이션** (Raised-cosine Edge Apodization)

## Repository Structure

```
simulators/
  bm_simulator.html    — 대화형 2D FFT 회절 시뮬레이터 (브라우저에서 바로 실행)
docs/
  physics.md           — 물리 수식 유도 및 이론적 배경
```

## Simulator

`simulators/bm_simulator.html`을 브라우저에서 열면 바로 사용 가능합니다.

### 기능

- **실시간 2D FFT** 기반 회절 패턴 시뮬레이션 (128x128 그리드)
- **7가지 파라미터** 슬라이더 조절:
  - 화소 피치 (30~100μm)
  - BM 폭 (5~30μm)
  - 파장 (380~780nm)
  - 위상 변조 깊이 (0~π)
  - 체크보드 ON/OFF
  - 랜덤 나노 디퓨저 강도 (0~1)
  - 에지 아포다이제이션 폭 (0~5μm)
- **6가지 프리셋**: Baseline / Horizontal / Checkerboard / Apodization / Random Diffuser / Full Combo
- **성능 지표**: H/V 회절 억제율(dB), DC 반사(dB), 총 산란 파워 비율
- 3패널 디스플레이: 실공간 픽셀 그리드 | 2D 회절 패턴 | 성능 메트릭

### 물리 모델

단위 셀 투과 함수:

```
t(x,y) = rect(x/a)·rect(y/b) · exp(i·φ(x,y))
```

체크보드 위상 변조:

```
φ(m,n) = π · [(m + n) mod 2]    (체크보드)
```

회절 패턴:

```
I(kx,ky) = |FT{t(x,y)}|²
```

완전 소멸 조건: `φ = π`, 체크보드 배열

## Key Results

| 모드 | H-회절 (dB) | V-회절 (dB) | DC (dB) |
|------|------------|------------|---------|
| Baseline | 0 | 0 | 0 |
| Horizontal | -15 | +2 | -5 |
| **Checkerboard** | **-22** | **-22** | **-18** |
| Full Combo | **-25** | **-25** | **-20** |

## Background

### 왜 원편광판을 없애고 싶은가?

- 원편광판은 OLED 자체 발광의 ~50%를 흡수 → 밝기/수명 손실
- 원편광판 제거 시 발광 효율 2배 → 같은 밝기에서 구동 전류 절반 → 수명 대폭 연장
- 하지만 제거 시 BM 주기 구조에 의한 **외부광 반사 회절** (무지개 패턴) 발생

### 제안하는 해법

BM 자체의 구조를 변형하여 회절 차수를 억제:

1. **체크보드 위상**: 인접 픽셀에 π 위상차 → 1차 회절 destructive interference
2. **랜덤 디퓨저**: 나노입자 분산으로 위상 랜덤화 → 회절 에너지를 넓게 분산
3. **에지 아포다이제이션**: BM 경계를 점진적으로 변화 → 고차 회절 억제

## Related Work

- SID Display Week 2026 Poster P-99: "Stress-Induced Birefringence Measurement for Crack Detection in Bessel Beam Laser-Drilled Glass" (동일 저자)
- 삼성디스플레이 OLED BM 구조 최적화 연구

## License

연구 목적 공개. 상업적 사용 시 문의 필요.

---
*© 2026 waterfirst*

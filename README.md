## About
- [A perceptual color space for image processing](https://bottosson.github.io/posts/oklab/)
- [Okhsv and Okhsl — Two new color spaces for color picking](https://bottosson.github.io/posts/colorpicker/)

## API
```php
// Oklab
OklabColor::oklab_to_linear_srgb(float $L, float $a, float $b): array // [float $r, float $g, float $b]
OklabColor::linear_srgb_to_oklab(float $r, float $g, float $b): array // [float $L, float $a, float $b]

// Okhsl
OklabColor::okhsl_to_srgb(float $h, float $s, float $l): array // [float $r, float $g, float $b]
OklabColor::srgb_to_okhsl(float $r, float $g, float $b): array // [float $h, float $s, float $l]

// Okhsv
OklabColor::okhsv_to_srgb(float $h, float $s, float $v): array // [float $r, float $g, float $b]
OklabColor::srgb_to_okhsv(float $r, float $g, float $b): array // [float $h, float $s, float $v]

// etc
OklabColor::gamut_clip_preserve_chroma(float $r, float $g, float $b): array
OklabColor::gamut_clip_project_to_0_5(float $r, float $g, float $b): array
OklabColor::gamut_clip_project_to_L_cusp(float $r, float $g, float $b): array
OklabColor::gamut_clip_adaptive_L0_0_5(float $r, float $g, float $b, float $alpha = 0.05): array
OklabColor::gamut_clip_adaptive_L0_L_cusp(float $r, float $g, float $b, float $alpha = 0.05): array
```
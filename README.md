## About
- [A perceptual color space for image processing](https://bottosson.github.io/posts/oklab/)
- [Okhsv and Okhsl — Two new color spaces for color picking](https://bottosson.github.io/posts/colorpicker/)

## Install
`composer require s13k/php_oklab_color`

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
```

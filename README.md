## About
- [A perceptual color space for image processing](https://bottosson.github.io/posts/oklab/)
- [Okhsv and Okhsl — Two new color spaces for color picking](https://bottosson.github.io/posts/colorpicker/)
- [Original C++ source (ok_color.h)](https://bottosson.github.io/misc/ok_color.h)

## Install
`composer require s13k/php_oklab_color`

## API
```php
// Oklab
Oklab\OkColor::oklab_to_linear_srgb(float $L, float $a, float $b): array // [float $r, float $g, float $b]
Oklab\OkColor::linear_srgb_to_oklab(float $r, float $g, float $b): array // [float $L, float $a, float $b]

// Okhsl
Oklab\OkColor::okhsl_to_srgb(float $h, float $s, float $l): array // [float $r, float $g, float $b]
Oklab\OkColor::srgb_to_okhsl(float $r, float $g, float $b): array // [float $h, float $s, float $l]

// Okhsv
Oklab\OkColor::okhsv_to_srgb(float $h, float $s, float $v): array // [float $r, float $g, float $b]
Oklab\OkColor::srgb_to_okhsv(float $r, float $g, float $b): array // [float $h, float $s, float $v]
```

## Example
```php
// Convert HSL (0–360°, 0–100%, 0–100%) to RGB hex
function okhsl_to_hex($h, $s, $l) {
    list($r, $g, $b) = Oklab\OkColor::okhsl_to_srgb($h % 360 / 360, $s / 100, $l / 100);
    return sprintf('#%02x%02x%02x', $r * 255, $g * 255, $b * 255);
}
echo okhsl_to_hex(137, 13, 37); #51594e
```

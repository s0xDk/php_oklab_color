<?php
namespace Oklab;

class OkColor
{

	private static function clamp(float $x, float $min, float $max): float
	{
		if ($x < $min)
			return $min;
		if ($x > $max)
			return $max;

		return $x;
	}

	private static function sgn(float $x): float
	{
		return (float)(0 < $x) - (float)($x < 0);
	}

	private static function srgb_transfer_function(float $a): float
	{
		return .0031308 >= $a ? 12.92 * $a : 1.055 * pow($a, .4166666666666667) - .055;
	}

	private static function srgb_transfer_function_inv(float $a): float
	{
		return .04045 < $a ? pow(($a + .055) / 1.055, 2.4) : $a / 12.92;
	}

	public static function linear_srgb_to_oklab(float $r, float $g, float $b): array
	{
		$l = 0.4122214708 * $r + 0.5363325363 * $g + 0.0514459929 * $b;
		$m = 0.2119034982 * $r + 0.6806995451 * $g + 0.1073969566 * $b;
		$s = 0.0883024619 * $r + 0.2817188376 * $g + 0.6299787005 * $b;

		$l_ = $l ** (1 / 3);
		$m_ = $m ** (1 / 3);
		$s_ = $s ** (1 / 3);

		return [
			0.2104542553 * $l_ + 0.7936177850 * $m_ - 0.0040720468 * $s_,
			1.9779984951 * $l_ - 2.4285922050 * $m_ + 0.4505937099 * $s_,
			0.0259040371 * $l_ + 0.7827717662 * $m_ - 0.8086757660 * $s_,
		];
	}

	public static function oklab_to_linear_srgb(float $L, float $a, float $b): array
	{
		$l_ = $L + 0.3963377774 * $a + 0.2158037573 * $b;
		$m_ = $L - 0.1055613458 * $a - 0.0638541728 * $b;
		$s_ = $L - 0.0894841775 * $a - 1.2914855480 * $b;

		$l = $l_ * $l_ * $l_;
		$m = $m_ * $m_ * $m_;
		$s = $s_ * $s_ * $s_;

		return [
			(+4.0767416621 * $l - 3.3077115913 * $m + 0.2309699292 * $s),
			(-1.2684380046 * $l + 2.6097574011 * $m - 0.3413193965 * $s),
			(-0.0041960863 * $l - 0.7034186147 * $m + 1.7076147010 * $s),
		];
	}

	// Finds the maximum saturation possible for a given hue that fits in sRGB
	// Saturation here is defined as S = C/L
	// a and b must be normalized so a^2 + b^2 == 1
	private static function compute_max_saturation(float $a, float $b): float
	{
		// Max saturation will be when one of r, g or b goes below zero.

		// Select different coefficients depending on which component goes below zero first

		if (-1.88170328 * $a - 0.80936493 * $b > 1)
		{
			// Red component
			$k0 = +1.19086277; $k1 = +1.76576728; $k2 = +0.59662641; $k3 = +0.75515197; $k4 = +0.56771245;
			$wl = +4.0767416621; $wm = -3.3077115913; $ws = +0.2309699292;
		}
		elseif (1.81444104 * $a - 1.19445276 * $b > 1)
		{
			// Green component
			$k0 = +0.73956515; $k1 = -0.45954404; $k2 = +0.08285427; $k3 = +0.12541070; $k4 = +0.14503204;
			$wl = -1.2684380046; $wm = +2.6097574011; $ws = -0.3413193965;
		}
		else
		{
			// Blue component
			$k0 = +1.35733652; $k1 = -0.00915799; $k2 = -1.15130210; $k3 = -0.50559606; $k4 = +0.00692167;
			$wl = -0.0041960863; $wm = -0.7034186147; $ws = +1.7076147010;
		}

		// Approximate max saturation using a polynomial:
		$S = $k0 + $k1 * $a + $k2 * $b + $k3 * $a * $a + $k4 * $a * $b;

		// Do one step Halley's method to get closer
		// this gives an error less than 10e6, except for some blue hues where the dS/dh is close to infinite
		// this should be sufficient for most applications, otherwise do two/three steps

		$k_l = +0.3963377774 * $a + 0.2158037573 * $b;
		$k_m = -0.1055613458 * $a - 0.0638541728 * $b;
		$k_s = -0.0894841775 * $a - 1.2914855480 * $b;

		{
			$l_ = 1 + $S * $k_l;
			$m_ = 1 + $S * $k_m;
			$s_ = 1 + $S * $k_s;

			$l = $l_ * $l_ * $l_;
			$m = $m_ * $m_ * $m_;
			$s = $s_ * $s_ * $s_;

			$l_dS = 3 * $k_l * $l_ * $l_;
			$m_dS = 3 * $k_m * $m_ * $m_;
			$s_dS = 3 * $k_s * $s_ * $s_;

			$l_dS2 = 6 * $k_l * $k_l * $l_;
			$m_dS2 = 6 * $k_m * $k_m * $m_;
			$s_dS2 = 6 * $k_s * $k_s * $s_;

			$f  = $wl * $l     + $wm * $m     + $ws * $s;
			$f1 = $wl * $l_dS  + $wm * $m_dS  + $ws * $s_dS;
			$f2 = $wl * $l_dS2 + $wm * $m_dS2 + $ws * $s_dS2;

			$S = $S - $f * $f1 / ($f1 * $f1 - 0.5 * $f * $f2);
		}

		return $S;

	}

	// finds L_cusp and C_cusp for a given hue
	// a and b must be normalized so a^2 + b^2 == 1
	private static function find_cusp(float $a, float $b): array
	{
		// First, find the maximum saturation (saturation S = C/L)
		$S_cusp = self::compute_max_saturation($a, $b);

		// Convert to linear sRGB to find the first point where at least one of r,g or b >= 1:
		$rgb_at_max = self::oklab_to_linear_srgb(1, $S_cusp * $a, $S_cusp * $b);
		$L_cusp = (1 / max(max($rgb_at_max[0], $rgb_at_max[1]), $rgb_at_max[2])) ** (1 / 3);
		$C_cusp = $L_cusp * $S_cusp;

		return [ $L_cusp, $C_cusp ];
	}

	// Finds intersection of the line defined by
	// L = L0 * (1 - t) + t * L1;
	// C = t * C1;
	// a and b must be normalized so a^2 + b^2 == 1
	private static function find_gamut_intersection(float $a, float $b, float $L1, float $C1, float $L0, array $cusp = null): float
	{
		// Find the intersection for upper and lower half seprately
		if ((($L1 - $L0) * $cusp[1] - ($cusp[0] - $L0) * $C1) <= 0)
		{
			// Lower half

			$t = $cusp[1] * $L0 / ($C1 * $cusp[0] + $cusp[1] * ($L0 - $L1));
		}
		else
		{
			// Upper half

			// First intersect with triangle
			$t = $cusp[1] * ($L0 - 1) / ($C1 * ($cusp[0] - 1) + $cusp[1] * ($L0 - $L1));

			// Then one step Halley's method
			{
				$dL = $L1 - $L0;
				$dC = $C1;

				$k_l = +0.3963377774 * $a + 0.2158037573 * $b;
				$k_m = -0.1055613458 * $a - 0.0638541728 * $b;
				$k_s = -0.0894841775 * $a - 1.2914855480 * $b;

				$l_dt = $dL + $dC * $k_l;
				$m_dt = $dL + $dC * $k_m;
				$s_dt = $dL + $dC * $k_s;

				// If higher accuracy is required, 2 or 3 iterations of the following block can be used:
				{
					$L = $L0 * (1 - $t) + $t * $L1;
					$C = $t * $C1;

					$l_ = $L + $C * $k_l;
					$m_ = $L + $C * $k_m;
					$s_ = $L + $C * $k_s;

					$l = $l_ * $l_ * $l_;
					$m = $m_ * $m_ * $m_;
					$s = $s_ * $s_ * $s_;

					$ldt = 3 * $l_dt * $l_ * $l_;
					$mdt = 3 * $m_dt * $m_ * $m_;
					$sdt = 3 * $s_dt * $s_ * $s_;

					$ldt2 = 6 * $l_dt * $l_dt * $l_;
					$mdt2 = 6 * $m_dt * $m_dt * $m_;
					$sdt2 = 6 * $s_dt * $s_dt * $s_;

					$r = 4.0767416621 * $l - 3.3077115913 * $m + 0.2309699292 * $s - 1;
					$r1 = 4.0767416621 * $ldt - 3.3077115913 * $mdt + 0.2309699292 * $sdt;
					$r2 = 4.0767416621 * $ldt2 - 3.3077115913 * $mdt2 + 0.2309699292 * $sdt2;

					$u_r = $r1 / ($r1 * $r1 - 0.5 * $r * $r2);
					$t_r = -$r * $u_r;

					$g = -1.2684380046 * $l + 2.6097574011 * $m - 0.3413193965 * $s - 1;
					$g1 = -1.2684380046 * $ldt + 2.6097574011 * $mdt - 0.3413193965 * $sdt;
					$g2 = -1.2684380046 * $ldt2 + 2.6097574011 * $mdt2 - 0.3413193965 * $sdt2;

					$u_g = $g1 / ($g1 * $g1 - 0.5 * $g * $g2);
					$t_g = -$g * $u_g;

					$b = -0.0041960863 * $l - 0.7034186147 * $m + 1.7076147010 * $s - 1;
					$b1 = -0.0041960863 * $ldt - 0.7034186147 * $mdt + 1.7076147010 * $sdt;
					$b2 = -0.0041960863 * $ldt2 - 0.7034186147 * $mdt2 + 1.7076147010  * $sdt2;

					$u_b = $b1 / ($b1 * $b1 - 0.5 * $b * $b2);
					$t_b = -$b * $u_b;

					$t_r = $u_r >= 0 ? $t_r : PHP_FLOAT_MAX;
					$t_g = $u_g >= 0 ? $t_g : PHP_FLOAT_MAX;
					$t_b = $u_b >= 0 ? $t_b : PHP_FLOAT_MAX;

					$t += min($t_r, min($t_g, $t_b));
				}
			}
		}

		return $t;
	}

	public static function gamut_clip_preserve_chroma(float $r, float $g, float $b): array
	{
		if ($r < 1 && $g < 1 && $b < 1 && $r > 0 && $g > 0 && $b > 0)
			return [$r, $g, $b];

		$lab = self::linear_srgb_to_oklab($r, $g, $b);

		$L = $lab[0];
		$eps = 0.00001;
		$C = max($eps, sqrt($lab[1] * $lab[1] + $lab[2] * $lab[2]));
		$a_ = $lab[1] / $C;
		$b_ = $lab[2] / $C;

		$L0 = self::clamp($L, 0, 1);

		$t = self::find_gamut_intersection($a_, $b_, $L, $C, $L0);
		$L_clipped = $L0 * (1 - $t) + $t * $L;
		$C_clipped = $t * $C;

		return self::oklab_to_linear_srgb($L_clipped, $C_clipped * $a_, $C_clipped * $b_);
	}

	public static function gamut_clip_project_to_0_5(float $r, float $g, float $b): array
	{
		if ($r < 1 && $g < 1 && $b < 1 && $r > 0 && $g > 0 && $b > 0)
			return [$r, $g, $b];

		$lab = self::linear_srgb_to_oklab($r, $g, $b);

		$L = $lab[0];
		$eps = 0.00001;
		$C = max($eps, sqrt($lab[1] * $lab[1] + $lab[2] * $lab[2]));
		$a_ = $lab[1] / $C;
		$b_ = $lab[2] / $C;

		$L0 = 0.5;

		$t = self::find_gamut_intersection($a_, $b_, $L, $C, $L0);
		$L_clipped = $L0 * (1 - $t) + $t * $L;
		$C_clipped = $t * $C;

		return self::oklab_to_linear_srgb($L_clipped, $C_clipped * $a_, $C_clipped * $b_);
	}

	public static function gamut_clip_project_to_L_cusp(float $r, float $g, float $b): array
	{
		if ($r < 1 && $g < 1 && $b < 1 && $r > 0 && $g > 0 && $b > 0)
			return [$r, $g, $b];

		$lab = self::linear_srgb_to_oklab($r, $g, $b);

		$L = $lab[0];
		$eps = 0.00001;
		$C = max($eps, sqrt($lab[1] * $lab[1] + $lab[2] * $lab[2]));
		$a_ = $lab[1] / $C;
		$b_ = $lab[2] / $C;

		// The cusp is computed here and in find_gamut_intersection, an optimized solution would only compute it once.
		$cusp = self::find_cusp($a_, $b_);

		$L0 = $cusp[0];

		$t = self::find_gamut_intersection($a_, $b_, $L, $C, $L0);
		$L_clipped = $L0 * (1 - $t) + $t * $L;
		$C_clipped = $t * $C;

		return self::oklab_to_linear_srgb($L_clipped, $C_clipped * $a_, $C_clipped * $b_);
	}

	public static function gamut_clip_adaptive_L0_0_5(float $r, float $g, float $b, float $alpha = 0.05): array
	{
		if ($r < 1 && $g < 1 && $b < 1 && $r > 0 && $g > 0 && $b > 0)
            return [$r, $g, $b];

		$lab = self::linear_srgb_to_oklab($r, $g, $b);

		$L = $lab[0];
		$eps = 0.00001;
		$C = max($eps, sqrt($lab[1] * $lab[1] + $lab[2] * $lab[2]));
		$a_ = $lab[1] / $C;
		$b_ = $lab[2] / $C;

		$Ld = $L - 0.5;
		$e1 = 0.5 + abs($Ld) + $alpha * $C;
		$L0 = 0.5 * (1 + self::sgn($Ld) * ($e1 - sqrt($e1 * $e1 - 2 * abs($Ld))));

		$t = self::find_gamut_intersection($a_, $b_, $L, $C, $L0);
		$L_clipped = $L0 * (1 - $t) + $t * $L;
		$C_clipped = $t * $C;

		return self::oklab_to_linear_srgb($L_clipped, $C_clipped * $a_, $C_clipped * $b_);
	}

	public static function gamut_clip_adaptive_L0_L_cusp(float $r, float $g, float $b, float $alpha = 0.05): array
	{
		if ($r < 1 && $g < 1 && $b < 1 && $r > 0 && $g > 0 && $b > 0)
            return [$r, $g, $b];

		$lab = self::linear_srgb_to_oklab($r, $g, $b);

		$L = $lab[0];
		$eps = 0.00001;
		$C = max($eps, sqrt($lab[1] * $lab[1] + $lab[2] * $lab[2]));
		$a_ = $lab[1] / $C;
		$b_ = $lab[2] / $C;

		// The cusp is computed here and in find_gamut_intersection, an optimized solution would only compute it once.
		$cusp = self::find_cusp($a_, $b_);

		$Ld = $L - $cusp[0];
		$k = 2 * ($Ld > 0 ? 1 - $cusp[0] : $cusp[0]);

		$e1 = 0.5 * $k + abs($Ld) + $alpha * $C / $k;
		$L0 = $cusp[0] + 0.5 * (self::sgn($Ld) * ($e1 - sqrt($e1 * $e1 - 2 * $k * abs($Ld))));

		$t = self::find_gamut_intersection($a_, $b_, $L, $C, $L0);
		$L_clipped = $L0 * (1 - $t) + $t * $L;
		$C_clipped = $t * $C;

		return self::oklab_to_linear_srgb($L_clipped, $C_clipped * $a_, $C_clipped * $b_);
	}

	private static function toe(float $x): float
	{
		$k_1 = 0.206;
		$k_2 = 0.03;
		$k_3 = (1 + $k_1) / (1 + $k_2);
		return 0.5 * ($k_3 * $x - $k_1 + sqrt(($k_3 * $x - $k_1) * ($k_3 * $x - $k_1) + 4 * $k_2 * $k_3 * $x));
	}

	private static function toe_inv(float $x): float
	{
		$k_1 = 0.206;
		$k_2 = 0.03;
		$k_3 = (1 + $k_1) / (1 + $k_2);
		return ($x * $x + $k_1 * $x) / ($k_3 * ($x + $k_2));
	}

	private static function to_ST(array $cusp): array
	{
		$L = $cusp[0];
		$C = $cusp[1];
		return [ $C / $L, $C / (1 - $L) ];
	}

	// Returns a smooth approximation of the location of the cusp
	// This polynomial was created by an optimization process
	// It has been designed so that S_mid < S_max and T_mid < T_max
	private static function get_ST_mid(float $a_, float $b_): array
	{
		$S = 0.11516993 + 1 / (
			+7.44778970 + 4.15901240 * $b_
			+ $a_ * (-2.19557347 + 1.75198401 * $b_
				+ $a_ * (-2.13704948 - 10.02301043 * $b_
					+ $a_ * (-4.24894561 + 5.38770819 * $b_ + 4.69891013 * $a_
						)))
			);

		$T = 0.11239642 + 1 / (
			+1.61320320 - 0.68124379 * $b_
			+ $a_ * (+0.40370612 + 0.90148123 * $b_
				+ $a_ * (-0.27087943 + 0.61223990 * $b_
					+ $a_ * (+0.00299215 - 0.45399568 * $b_ - 0.14661872 * $a_
						)))
			);

		return [ $S, $T ];
	}

	private static function get_Cs(float $L, float $a_, float $b_): array
	{
		$cusp = self::find_cusp($a_, $b_);

		$C_max = self::find_gamut_intersection($a_, $b_, $L, 1, $L, $cusp);
		$ST_max = self::to_ST($cusp);

		// Scale factor to compensate for the curved part of gamut shape:
		$k = $C_max / min(($L * $ST_max[0]), (1 - $L) * $ST_max[1]);

		$C_mid = 0;
		{
			$ST_mid = self::get_ST_mid($a_, $b_);

			// Use a soft minimum function, instead of a sharp triangle shape to get a smooth value for chroma.
			$C_a = $L * $ST_mid[0];
			$C_b = (1 - $L) * $ST_mid[1];
			$C_mid = 0.9 * $k * sqrt(sqrt(1 / (1 / ($C_a * $C_a * $C_a * $C_a) + 1 / ($C_b * $C_b * $C_b * $C_b))));
		}

		$C_0 = 0;
		{
			// for C_0, the shape is independent of hue, so ST are constant. Values picked to roughly be the average values of ST.
			$C_a = $L * 0.4;
			$C_b = (1 - $L) * 0.8;

			// Use a soft minimum function, instead of a sharp triangle shape to get a smooth value for chroma.
			$C_0 = sqrt(1 / (1 / ($C_a * $C_a) + 1 / ($C_b * $C_b)));
		}

		return [ $C_0, $C_mid, $C_max ];
	}

	public static function okhsl_to_srgb(float $h, float $s, float $l): array
	{
		if ($l == 1)
		{
			return [1, 1, 1];
		}
		elseif ($l == 0)
		{
			return [0, 0, 0];
		}

		$a_ = cos(2 * M_PI * $h);
		$b_ = sin(2 * M_PI * $h);
		$L = self::toe_inv($l);

		$Cs = self::get_Cs($L, $a_, $b_);
		$C_0 = $Cs[0];
		$C_mid = $Cs[1];
		$C_max = $Cs[2];

		$mid = 0.8;
		$mid_inv = 1.25;

		if ($s < $mid)
		{
			$t = $mid_inv * $s;

			$k_1 = $mid * $C_0;
			$k_2 = (1 - $k_1 / $C_mid);

			$C = $t * $k_1 / (1 - $k_2 * $t);
		}
		else
		{
			$t = 5 * ($s - $mid);

			$k_0 = $C_mid;
			$k_1 = (1 - $mid) * $C_mid * $C_mid * $mid_inv * $mid_inv / $C_0;
			$k_2 = (1 - ($k_1) / ($C_max - $C_mid));

			$C = $k_0 + $t * $k_1 / (1 - $k_2 * $t);
		}

		$rgb = self::oklab_to_linear_srgb($L, $C * $a_, $C * $b_);
		return [
			self::srgb_transfer_function($rgb[0]),
			self::srgb_transfer_function($rgb[1]),
			self::srgb_transfer_function($rgb[2]),
		];
	}

	public static function srgb_to_okhsl(float $r, float $g, float $b): array
	{
		$lab = self::linear_srgb_to_oklab(
			self::srgb_transfer_function_inv($r),
			self::srgb_transfer_function_inv($g),
			self::srgb_transfer_function_inv($b)
		);

		$C = sqrt($lab[1] * $lab[1] + $lab[2] * $lab[2]);
		$a_ = $lab[1] / $C;
		$b_ = $lab[2] / $C;

		$L = $lab[0];
		$h = 0.5 + 0.5 * atan2(-$lab[2], -$lab[1]) / M_PI;

		$cs = self::get_Cs($L, $a_, $b_);
		$C_0 = $cs[0];
		$C_mid = $cs[1];
		$C_max = $cs[2];

		// Inverse of the interpolation in okhsl_to_srgb:

		$mid = 0.8;
		$mid_inv = 1.25;

		if ($C < $C_mid)
		{
			$k_1 = $mid * $C_0;
			$k_2 = (1 - $k_1 / $C_mid);

			$t = $C / ($k_1 + $k_2 * $C);
			$s = $t * $mid;
		}
		else
		{
			$k_0 = $C_mid;
			$k_1 = (1 - $mid) * $C_mid * $C_mid * $mid_inv * $mid_inv / $C_0;
			$k_2 = (1 - ($k_1) / ($C_max - $C_mid));

			$t = ($C - $k_0) / ($k_1 + $k_2 * ($C - $k_0));
			$s = $mid + (1 - $mid) * $t;
		}

		$l = self::toe($L);
		return [ $h, $s, $l ];
	}

	public static function okhsv_to_srgb(float $h, float $s, float $v): array
	{
		$a_ = cos(2 * M_PI * $h);
		$b_ = sin(2 * M_PI * $h);

		$cusp = self::find_cusp($a_, $b_);
		$ST_max = self::to_ST($cusp);
		$S_max = $ST_max[0];
		$T_max = $ST_max[1];
		$S_0 = 0.5;
		$k = 1 - $S_0 / $S_max;

		// first we compute L and V as if the gamut is a perfect triangle:

		// L, C when v==1:
		$L_v = 1      - $s * $S_0 / ($S_0 + $T_max - $T_max * $k * $s);
		$C_v = $s * $T_max * $S_0 / ($S_0 + $T_max - $T_max * $k * $s);

		$L = $v * $L_v;
		$C = $v * $C_v;

		// then we compensate for both toe and the curved top part of the triangle:
		$L_vt = self::toe_inv($L_v);
		$C_vt = $C_v * $L_vt / $L_v;

		$L_new = self::toe_inv($L);
		$C = $C * $L_new / $L;
		$L = $L_new;

		$rgb_scale = self::oklab_to_linear_srgb($L_vt, $a_ * $C_vt, $b_ * $C_vt);
		$scale_L = (1 / max(max($rgb_scale[0], $rgb_scale[1]), max($rgb_scale[2], 0))) ** (1/3);

		$L = $L * $scale_L;
		$C = $C * $scale_L;

		$rgb = self::oklab_to_linear_srgb($L, $C * $a_, $C * $b_);
		return [
			self::srgb_transfer_function($rgb[0]),
			self::srgb_transfer_function($rgb[1]),
			self::srgb_transfer_function($rgb[2]),
		];
	}

	public static function srgb_to_okhsv(float $r, float $g, float $b): array
	{
		$lab = self::linear_srgb_to_oklab(
			self::srgb_transfer_function_inv($r),
			self::srgb_transfer_function_inv($g),
			self::srgb_transfer_function_inv($b)
		);

		$C = sqrt($lab[1] * $lab[1] + $lab[2] * $lab[2]);
		$a_ = $lab[1] / $C;
		$b_ = $lab[2] / $C;

		$L = $lab[0];
		$h = 0.5 + 0.5 * atan2(-$lab[2], -$lab[1]) / M_PI;

		$cusp = self::find_cusp($a_, $b_);
		$ST_max = self::to_ST($cusp);
		$S_max = $ST_max[0];
		$T_max = $ST_max[1];
		$S_0 = 0.5;
		$k = 1 - $S_0 / $S_max;

		// first we find L_v, C_v, L_vt and C_vt

		$t = $T_max / ($C + $L * $T_max);
		$L_v = $t * $L;
		$C_v = $t * $C;

		$L_vt = self::toe_inv($L_v);
		$C_vt = $C_v * $L_vt / $L_v;

		// we can then use these to invert the step that compensates for the toe and the curved top part of the triangle:
		$rgb_scale = self::oklab_to_linear_srgb($L_vt, $a_ * $C_vt, $b_ * $C_vt);
		$scale_L = (1 / max(max($rgb_scale[0], $rgb_scale[1]), max($rgb_scale[2], 0))) ** (1/3);

		$L = $L / $scale_L;
		$C = $C / $scale_L;

		$C = $C * self::toe($L) / $L;
		$L = self::toe($L);

		// we can now compute v and s:

		$v = $L / $L_v;
		$s = ($S_0 + $T_max) * $C_v / (($T_max * $S_0) + $T_max * $k * $C_v);

		return [ $h, $s, $v ];
	}

}

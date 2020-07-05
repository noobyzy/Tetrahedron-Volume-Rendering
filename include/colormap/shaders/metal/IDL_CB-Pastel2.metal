#include <metal_stdlib>
using namespace metal;

namespace colormap {
namespace IDL {
namespace CB_Pastel2 {

float colormap_red(float x) {
	if (x < 0.1414470427532423) {
		return 5.23716927453769E+02 * x + 1.79102418207681E+02;
	} else if (x < 0.2832126252873305) {
		return -3.55011583011593E+02 * x + 3.03395967395968E+02;
	} else if (x < 0.4293789173286286) {
		return 2.81737389211071E+02 * x + 1.23060619323778E+02;
	} else if (x < 0.5703484841123749) {
		return -9.85406162465110E+01 * x + 2.86343977591045E+02;
	} else if (x < 0.7170614267751989) {
		return 1.69092460881909E+02 * x + 1.33699857752520E+02;
	} else if (x < 0.859829619768543) {
		return -9.94710581026121E+01 * x + 3.26276397855329E+02;
	} else {
		return -2.57056149732620E+02 * x + 4.61772727272750E+02;
	}
}

float colormap_green(float x) {
	if (x < 0.1411063922737659) {
		return -1.49220483641537E+02 * x + 2.26007112375533E+02;
	} else if (x < 0.2816283290590322) {
		return 5.77629343629288E+01 * x + 1.96800429000430E+02;
	} else if (x < 0.4291887492428612) {
		return -7.38876244665610E+01 * x + 2.33876955903267E+02;
	} else if (x < 0.571830104540257) {
		return 3.01873399715509E+02 * x + 7.26045519203479E+01;
	} else if (x < 0.7190262682310248) {
		return -2.25206477732972E+01 * x + 2.58102834008109E+02;
	} else if (x < 0.8491803538380496) {
		return -1.16468292682893E+02 * x + 3.25653658536549E+02;
	} else {
		return -1.44447728516695E+02 * x + 3.49413245758086E+02;
	}
}

float colormap_blue(float x) {
	if (x < 0.1425168965591466) {
		return -2.26903238866400E+02 * x + 2.04742307692308E+02;
	} else if (x < 0.2851292529683606) {
		return 4.18120091673021E+02 * x + 1.12815584415585E+02;
	} else if (x < 0.4319360871262316) {
		return -3.09335813546247E+01 * x + 2.40853922748656E+02;
	} else if (x < 0.7146533590447866) {
		return -1.88956299440485E+02 * x + 3.09109637275714E+02;
	} else if (x < 0.8619542566532371) {
		return 2.06196082722327E+02 * x + 2.67126600285119E+01;
	} else {
		return -6.48097784562050E+00 * x + 2.10030557677552E+02;
	}
}

float4 colormap(float x) {
	float r = clamp(colormap_red(x) / 255.0, 0.0, 1.0);
	float g = clamp(colormap_green(x) / 255.0, 0.0, 1.0);
	float b = clamp(colormap_blue(x) / 255.0, 0.0, 1.0);
	return float4(r, g, b, 1.0);
}

} // namespace CB_Pastel2
} // namespace IDL
} // namespace colormap

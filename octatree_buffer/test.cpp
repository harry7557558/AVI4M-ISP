// test octatree.h + sample rendering

#pragma once

#define _MULTITHREADING 1

#pragma warning(disable: 4244)
#pragma warning(disable: 4305)

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <chrono>
#include <thread>
#include <string>

#include <Windows.h>
#include <windowsx.h>
#include <tchar.h>


#define WIN_NAME "Shadertoy Template"
#define WinW_Padding 100
#define WinH_Padding 100
#define WinW_Default 400
#define WinH_Default 400
#define WinW_Min 100
#define WinH_Min 100
#define WinW_Max 800
#define WinH_Max 800
// for my version of Windows
#define WinW_Offset 16
#define WinH_Offset 39

void Init(char* argv[]);
void render();
void WindowResize(int _oldW, int _oldH, int _W, int _H);
void WindowClose();
void MouseMove(int _X, int _Y);
void MouseWheel(int _DELTA);
void MouseDownL(int _X, int _Y);
void MouseUpL(int _X, int _Y);
void MouseDownR(int _X, int _Y);
void MouseUpR(int _X, int _Y);
void KeyDown(WPARAM _KEY);
void KeyUp(WPARAM _KEY);

HWND _HWND; int _WIN_W, _WIN_H;
HBITMAP _HIMG; COLORREF *_WINIMG;
#define Canvas(x,y) _WINIMG[(y)*_WIN_W+(x)]
#define setColor(x,y,col) do{if((x)>=0&&(x)<_WIN_W&&(y)>=0&&(y)<_WIN_H)Canvas(x,y)=col;}while(0)

auto _TimeStart = std::chrono::high_resolution_clock::now();
bool _RenderNeeded = true;


// Win32 Entry

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam) {
#define _RDBK { if (!_RenderNeeded) break; HDC hdc = GetDC(hWnd), HImgMem = CreateCompatibleDC(hdc); HBITMAP hbmOld = (HBITMAP)SelectObject(HImgMem, _HIMG); render(); BitBlt(hdc, 0, 0, _WIN_W, _WIN_H, HImgMem, 0, 0, SRCCOPY); SelectObject(HImgMem, hbmOld), DeleteDC(HImgMem), DeleteDC(hdc); _RenderNeeded = false; break; }
	switch (message) {
	case WM_NULL: { _RDBK }
	case WM_CREATE: { break; } case WM_CLOSE: { WindowClose(); DestroyWindow(hWnd); return 0; } case WM_DESTROY: { PostQuitMessage(0); return 0; }
	case WM_MOVE:; case WM_SIZE: {
		RECT Client; GetClientRect(hWnd, &Client);
		WindowResize(_WIN_W, _WIN_H, Client.right, Client.bottom);
		_WIN_W = Client.right, _WIN_H = Client.bottom;
		BITMAPINFO bmi; bmi.bmiHeader.biSize = sizeof(BITMAPINFO), bmi.bmiHeader.biWidth = Client.right, bmi.bmiHeader.biHeight = Client.bottom, bmi.bmiHeader.biPlanes = 1, bmi.bmiHeader.biBitCount = 32; bmi.bmiHeader.biCompression = BI_RGB, bmi.bmiHeader.biSizeImage = 0, bmi.bmiHeader.biXPelsPerMeter = bmi.bmiHeader.biYPelsPerMeter = 0, bmi.bmiHeader.biClrUsed = bmi.bmiHeader.biClrImportant = 0; bmi.bmiColors[0].rgbBlue = bmi.bmiColors[0].rgbGreen = bmi.bmiColors[0].rgbRed = bmi.bmiColors[0].rgbReserved = 0;
		if (_HIMG != NULL) DeleteObject(_HIMG);
		HDC hdc = GetDC(hWnd);
		_HIMG = CreateDIBSection(hdc, &bmi, DIB_RGB_COLORS, (void**)&_WINIMG, NULL, 0);
		DeleteDC(hdc);
		_RDBK }
	case WM_GETMINMAXINFO: { LPMINMAXINFO lpMMI = (LPMINMAXINFO)lParam; lpMMI->ptMinTrackSize.x = WinW_Min + WinW_Offset, lpMMI->ptMinTrackSize.y = WinH_Min + WinH_Offset, lpMMI->ptMaxTrackSize.x = WinW_Max + WinW_Offset, lpMMI->ptMaxTrackSize.y = WinH_Max + WinH_Offset; break; }
	case WM_PAINT: { PAINTSTRUCT ps; HDC hdc = BeginPaint(hWnd, &ps), HMem = CreateCompatibleDC(hdc); HBITMAP hbmOld = (HBITMAP)SelectObject(HMem, _HIMG); BitBlt(hdc, 0, 0, _WIN_W, _WIN_H, HMem, 0, 0, SRCCOPY); SelectObject(HMem, hbmOld); EndPaint(hWnd, &ps); DeleteDC(HMem), DeleteDC(hdc); break; }
#define _USER_FUNC_PARAMS GET_X_LPARAM(lParam), _WIN_H - 1 - GET_Y_LPARAM(lParam)
	case WM_MOUSEMOVE: { MouseMove(_USER_FUNC_PARAMS); _RDBK } case WM_MOUSEWHEEL: { MouseWheel(GET_WHEEL_DELTA_WPARAM(wParam)); _RDBK }
	case WM_LBUTTONDOWN: { SetCapture(hWnd); MouseDownL(_USER_FUNC_PARAMS); _RDBK } case WM_LBUTTONUP: { ReleaseCapture(); MouseUpL(_USER_FUNC_PARAMS); _RDBK }
	case WM_RBUTTONDOWN: { MouseDownR(_USER_FUNC_PARAMS); _RDBK } case WM_RBUTTONUP: { MouseUpR(_USER_FUNC_PARAMS); _RDBK }
	case WM_SYSKEYDOWN:; case WM_KEYDOWN: { if (wParam >= 0x08) KeyDown(wParam); _RDBK } case WM_SYSKEYUP:; case WM_KEYUP: { if (wParam >= 0x08) KeyUp(wParam); _RDBK }
	} return DefWindowProc(hWnd, message, wParam, lParam);
}
int main(int argc, char* argv[]) {
	Init(argv);
	HINSTANCE hInstance = NULL; int nCmdShow = SW_RESTORE;
	WNDCLASSEX wc;
	wc.cbSize = sizeof(WNDCLASSEX), wc.style = 0, wc.lpfnWndProc = WndProc, wc.cbClsExtra = wc.cbWndExtra = 0, wc.hInstance = hInstance;
	wc.hIcon = wc.hIconSm = 0, wc.hCursor = LoadCursor(NULL, IDC_ARROW), wc.hbrBackground = CreateSolidBrush(RGB(0, 0, 0)), wc.lpszMenuName = NULL, wc.lpszClassName = _T(WIN_NAME);
	if (!RegisterClassEx(&wc)) return -1;
	_HWND = CreateWindow(_T(WIN_NAME), _T(WIN_NAME), WS_OVERLAPPEDWINDOW, WinW_Padding, WinH_Padding, WinW_Default + WinW_Offset, WinH_Default + WinH_Offset, NULL, NULL, hInstance, NULL);
	ShowWindow(_HWND, nCmdShow); UpdateWindow(_HWND);
	MSG message; while (GetMessage(&message, 0, 0, 0)) {
		TranslateMessage(&message); DispatchMessage(&message);
	} return (int)message.wParam;
}



#include "../octatree/ply_writer.h"

#include "octatree.h"

#include "../octatree/trigs2mesh.h"


#define NAMESPACE_GLSL_BEGIN namespace GLSL {
#define NAMESPACE_GLSL_END }

// for attributing arrays to fit GLSL
#define out
#define inout

NAMESPACE_GLSL_BEGIN

float iRx = 0.3 + 1e-4;
float iRz = 0.5 + 1e-4;
float iSc = 0.5;
vec3 iResolution = vec3(0, 0, 1);
vec4 iMouse = vec4(0, 0, 0, 0);

vec4 *FrameBuffer = nullptr;

#define ZERO 0

#define P0 vec3(-2.5, -2.5, -2.0) /* min coordinates of grid */
#define P1 vec3(2.5, 2.5, 2.0) /* max coordinates of grid */
#define GRID_DIF ivec3(1) /* initial grid size, at least one odd component */
#define PLOT_DEPTH 9 /* depth of the tree */
#define GRID_SIZE (GRID_DIF*(1<<PLOT_DEPTH))
#define EDGE_ROUNDING 255 /* divide edge into # intervals and round to integer coordinate */
#define MESH_SIZE (GRID_SIZE*EDGE_ROUNDING)

#define GRID_EXPAND 4 /* pre-sample this number of layers in tree */
#define SEARCH_DIF_EXP (GRID_DIF*(1<<GRID_EXPAND))
#define PLOT_DEPTH_EXP (PLOT_DEPTH-GRID_EXPAND)


#if 1
namespace GLSL_INNER {
#define iTime 0.0
#define iFrame 0
#include "../.glsl.cpp"
};
vec4 map(vec3 p, bool col_required) {
	return GLSL_INNER::map(p, col_required);
}
#else
vec4 map(vec3 p, bool col_required) {
	float d = length(p) - 1.0 + 0.2*sin(10.0*p.x)*sin(10.0*p.y)*sin(10.0*p.z);
	//float d = min(length(vec2(length(p.xy()) - 1.0f, p.z)) - 0.5f, length(p) - 0.1f);
	//float d = sin(1.57*p.x)*sin(1.57*p.y)*sin(1.57*p.z) - 0.1;
	//float d = p.x + p.y - p.x*p.y;
	//float d = log(0.1*(exp(p.x) + exp(p.y) + exp(p.z)));
	//float d = max(abs(length(vec2(length(p.xy()) - 1.0, p.z)) - 0.5) - 0.05, p.z);
	return vec4(0.5 + 0.5*p, d);
}
#endif
float sdf(vec3 p) { return map(p, false).w; }
vec3 colorf(vec3 p) { return map(p, true).xyz(); }

//===========================

std::vector<uint8_t> treeSampler;
int getUint8(int i) {
	if (i >= (int)treeSampler.size()) {
		printf("%d %d\n", i, (int)treeSampler.size());
		return 0;
	}
	return treeSampler[i];
}
int getUint16(int i) {
	return getUint8(i) + getUint8(i + 1) * 256;
}
int getUint32(int i) {
	return getUint8(i) + 256 * (getUint8(i + 1) + 256 * (getUint8(i + 2) + 256 * getUint8(i + 3)));
}
ivec3 getUvec3(int i) {
	int x = getUint8(i);
	int y = getUint8(i + 1);
	int z = getUint8(i + 2);
	return ivec3(x, y, z);
}


#include "intersector_2.h"


#if 0
// fast preview
vec3 mainRender(vec3 ro, vec3 rd) {
	float t;
	vec3 n, col = vec3(0.0);
	if (intersectObject(ro, rd, t, 1e6, n, col)) {
		//return col;
		return col * abs(dot(normalize(n), rd));
	}
	return col;
}
#else
// ray tracing
vec3 mainRender(vec3 ro, vec3 rd) {
	const int MAT_BACKGROUND = 0;
	const int MAT_PLANE = 1;
	const int MAT_OBJECT = 2;

	vec3 m_col = vec3(1.0), t_col = vec3(0.0);
	bool inside_object = false;

	for (int iter = ZERO; iter < 64; iter++) {
		rd = normalize(rd);
		ro += 1e-4 * rd;
		float t, min_t = 1e12;
		vec3 n, min_n;
		int material = MAT_BACKGROUND;
		vec3 col = vec3(0.0);

		// plane
		t = -(ro.z - (P0.z - 0.01)) / rd.z;
		if (t > 0.0) {
			min_t = t, min_n = vec3(0, 0, 1);
			material = MAT_PLANE;
		}

		// object
		t = min_t;
		if (intersectObject(ro, rd, t, min_t, n, col) && t < min_t) {
			min_t = t, min_n = normalize(n);
			material = MAT_OBJECT;
		}

		// update ray
		if (material == MAT_BACKGROUND) {
			col = vec3(1.0)*max(rd.z, 0.0);
			return m_col * col + t_col;
		}
		min_n = dot(rd, min_n) < 0.0 ? min_n : -min_n;  // ray hits into the surface
		ro += rd * min_t;
		if (material == MAT_PLANE) {
			m_col *= vec3(0.8, 0.9, 1.0);
			rd = rd - 2.0*dot(rd, min_n)*min_n;
		}
		if (material == MAT_OBJECT) {
			m_col *= col;
			rd = rd - 2.0*dot(rd, min_n)*min_n;
		}
		if (inside_object) return 1e12*vec3(1, -1, -1);  // red warning
	}
	return m_col + t_col;
}
#endif

void mainImage(vec4 &fragColor, vec2 fragCoord) {
	//fragCoord = vec2(200.5, 234.5);
	const vec3 CENTER = vec3(0, 0, 0);
	const float DIST = 20.0f;  // larger = smaller
	const float VIEW_FIELD = 0.4f;  // larger = larger + more perspective

	// camera
	vec3 w = vec3(cos(iRx)*vec2(cos(iRz), sin(iRz)), sin(iRx));
	vec3 u = vec3(-sin(iRz), cos(iRz), 0);
	vec3 v = cross(w, u);
	vec3 ro = DIST * w + CENTER;
	vec2 uv = iSc * (2.0f*(fragCoord.xy() - 0.5f) / iResolution.xy() - 1.0f);
	vec2 sc = iResolution.xy() / min(iResolution.x, iResolution.y);
	vec3 rd = mat3(u, v, -w)*vec3(VIEW_FIELD*uv*sc, 1.0f);
	rd = normalize(rd);

	// calculate pixel color
	vec3 col = mainRender(ro, rd);
	fragColor = vec4(col, 1.0);
}

NAMESPACE_GLSL_END



void batchRender(void(*task)(int, int, int, bool*), int Max) {
#if _MULTITHREADING
	const int MAX_THREADS = std::thread::hardware_concurrency();
	bool* fn = new bool[MAX_THREADS];
	std::thread** T = new std::thread*[MAX_THREADS];
	for (int i = 0; i < MAX_THREADS; i++) {
		fn[i] = false;
		T[i] = new std::thread(task, i, Max, MAX_THREADS, &fn[i]);
	}
	int count; do {
		count = 0;
		for (int i = 0; i < MAX_THREADS; i++) count += fn[i];
	} while (count < MAX_THREADS);
	//for (int i = 0; i < MAX_THREADS; i++) delete T[i];
	delete fn; delete T;
#else
	task(0, Max, 1, NULL);
#endif
}

void render() {

	// initialize window
	for (int i = 0, l = _WIN_W * _WIN_H; i < l; i++) _WINIMG[i] = 0;

	// update buffer size
	int old_w = (int)GLSL::iResolution.x, old_h = (int)GLSL::iResolution.y;
	if (old_w != _WIN_W || old_h != _WIN_H) {
		if (GLSL::FrameBuffer) delete GLSL::FrameBuffer;
		GLSL::FrameBuffer = new vec4[_WIN_W*_WIN_H];
	}

	// set uniform variables
	GLSL::iResolution.x = (float)_WIN_W;
	GLSL::iResolution.y = (float)_WIN_H;

	// render frame
	auto time_0 = std::chrono::high_resolution_clock::now();
	batchRender([](int beg, int end, int step, bool* sig) {
		int WIN_SIZE = _WIN_W * _WIN_H;
		for (int k = beg; k < end; k += step) {
			int j = k / _WIN_W, i = k % _WIN_W;
			vec4 gl_FragCoord = vec4(vec2(i, j) + 0.5f, 0.0f, 0.0f);
			GLSL::mainImage(GLSL::FrameBuffer[k], gl_FragCoord.xy());
		}
		if (sig) *sig = true;
	}, _WIN_W*_WIN_H);
	float time_elapsed = std::chrono::duration<float>(std::chrono::high_resolution_clock::now() - time_0).count();

	// display image + statistics
	vec3 avrcol = vec3(0.0), varcol = vec3(0.0);
	vec3 mincol = vec3(1e12), maxcol = vec3(-1e12);
	for (int j = 0; j < _WIN_H; j++) {
		for (int i = 0; i < _WIN_W; i++) {
			vec3 rgb = GLSL::FrameBuffer[j*_WIN_W + i].xyz();
			avrcol += rgb, varcol += rgb * rgb, mincol = min(mincol, rgb), maxcol = max(maxcol, rgb);
			ivec3 irgb = ivec3(255.0*saturate(rgb) + 0.5);
			COLORREF c = (COLORREF(irgb.x) << 16) | (COLORREF(irgb.y) << 8) | (COLORREF(irgb.z));
			_WINIMG[j*_WIN_W + i] = c;
		}
	}
	int n = _WIN_W * _WIN_H;
	avrcol = avrcol / n;
	varcol = sqrt((varcol - n * avrcol * avrcol) / (n - 1));
	printf("mean=(%.3g,%.3g,%.3g), var=(%.3g,%.3g,%.3g), min=(%.3g,%.3g,%.3g), max=(%.3g,%.3g,%.3g)\n",
		avrcol.x, avrcol.y, avrcol.z, varcol.x, varcol.y, varcol.z, mincol.x, mincol.y, mincol.z, maxcol.x, maxcol.y, maxcol.z);

	// the actual fps is less because of display time
	char text[1024];
	sprintf(text, "%.1fms (%.1f fps) [%dx%d]", 1000.0f * time_elapsed, 1.0f / time_elapsed, _WIN_W, _WIN_H);
	SetWindowTextA(_HWND, text);
}



void exportModel(const char* filepath) {
	// compute parameters
	vec3 cell_size = (P1 - P0) / (vec3(SEARCH_DIF_EXP) * exp2(PLOT_DEPTH_EXP));
	float epsilon = 0.001f * std::min({ cell_size.x, cell_size.y, cell_size.z });
	printf("%d %d %d  %d  epsilon=%.2g\n", SEARCH_DIF_EXP.x, SEARCH_DIF_EXP.y, SEARCH_DIF_EXP.z, 1 << PLOT_DEPTH_EXP, epsilon);

	// marching cube
	int eval_count = 0;
	auto triangulate_start = std::chrono::high_resolution_clock::now();
	std::vector<triangle_3d> trigs = ScalarFieldTriangulator_octatree::octatree(
		[&](vec3 p) { eval_count++; return GLSL::sdf(p); },
		P0, P1,
		SEARCH_DIF_EXP, PLOT_DEPTH_EXP
	);
	float triangulate_time = std::chrono::duration<float>(std::chrono::high_resolution_clock::now() - triangulate_start).count();
	printf("%.1fms, %d evaluations\n", 1000.0f*triangulate_time, eval_count);
	printf("%d triangles => ", (int)trigs.size());

	// triangles to mesh
	std::vector<vec3> vertices;
	std::vector<ivec3> faces;
	TrigsToMesh(trigs, epsilon, vertices, faces);
	printf("%d vertices, %d faces\n", (int)vertices.size(), (int)faces.size());

	// color mesh
	auto coloring_start = std::chrono::high_resolution_clock::now();
	std::vector<vec3> colors;
	colors.resize((int)vertices.size());
	for (int i = 0; i < (int)vertices.size(); i++) {
		colors[i] = GLSL::colorf(vertices[i]);
	}
	float coloring_time = std::chrono::duration<float>(std::chrono::high_resolution_clock::now() - coloring_start).count();
	printf("%.1fms coloring\n", 1000.0f*coloring_time);

	// output
	WritePLY(filepath,
		&vertices[0], (int)vertices.size(),
		&faces[0], (int)faces.size(),
		&colors[0]
	);
}

void exportTree() {
	GLSL::treeSampler = ScalarFieldTriangulator_octatree::octatree_buffer(
		GLSL::sdf,
		P0, P1,
		SEARCH_DIF_EXP, PLOT_DEPTH_EXP, EDGE_ROUNDING,
		GLSL::colorf
	);
	printf("Tree sampler: %.2f MB\n", GLSL::treeSampler.size() / exp2(20));

	FILE* fp = fopen("D:\\.bin", "wb");
	fwrite(&GLSL::treeSampler[0], 1, GLSL::treeSampler.size(), fp);
	fclose(fp);
}

void Init(char* argv[]) {
#if 1
	exportModel(argv[1]); exit(0);
	//exportTree();
#else
	FILE* fp = fopen("D:\\.bin", "r");
	fseek(fp, 0, SEEK_END);
	int size = ftell(fp);
	rewind(fp);
	GLSL::treeSampler.resize(size);
	fread(&GLSL::treeSampler[0], 1, size, fp);
	fclose(fp);
#endif
}

void WindowResize(int _oldW, int _oldH, int _W, int _H) {

	_RenderNeeded = true;
}
void WindowClose() {
}

void MouseWheel(int _DELTA) {
	GLSL::iSc *= exp(-0.0005*_DELTA);
	_RenderNeeded = true;
}
void MouseDownL(int _X, int _Y) {
	GLSL::iMouse.zw() = vec2(_X, _Y);
	_RenderNeeded = true;
}
void MouseMove(int _X, int _Y) {
	static int old_x = -1;
	static int old_y = -1;
	if (GLSL::iMouse.z > 0.0f) {
		GLSL::iMouse.w = -abs(GLSL::iMouse.w);
		GLSL::iMouse.xy() = vec2(_X, _Y);
		if (old_x >= 0 && old_y >= 0) {
			int dx = _X - old_x;
			int dy = _Y - old_y;
			GLSL::iRz -= 4.0 * dx / GLSL::iResolution.x;
			GLSL::iRx -= 4.0 * dy / GLSL::iResolution.y;
		}
		old_x = _X, old_y = _Y;
	}
	else old_x = old_y = -1;
	_RenderNeeded = true;
}
void MouseUpL(int _X, int _Y) {
	printf("Mouse click at (%d, %d)\n", _X, _Y);
	GLSL::iMouse.zw() = -abs(GLSL::iMouse.zw());
	_RenderNeeded = true;
}
void MouseDownR(int _X, int _Y) {
	_RenderNeeded = true;
}
void MouseUpR(int _X, int _Y) {
	bool topmost = GetWindowLong(_HWND, GWL_EXSTYLE) & WS_EX_TOPMOST;
	SetWindowPos(_HWND, topmost ? HWND_NOTOPMOST : HWND_TOPMOST, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE);
	_RenderNeeded = true;
}
void KeyDown(WPARAM _KEY) {
}
void KeyUp(WPARAM _KEY) {
	_RenderNeeded = true;
}

// Place the following functions inside `#ifndef __CPLUSPLUS` in .glsl.cpp:

// bool intersectGlass(vec3 ro, vec3 rd, float &t, float t1, vec3 &n);
// bool intersectContent(vec3 ro, vec3 rd, float &t, float t1, vec3 &n, vec3 &col);


#define _MULTITHREADING 1


// data type conversion warning
#pragma warning(disable: 4244)
// data type truncation warning
#pragma warning(disable: 4305)

#include <cmath>
#include <stdio.h>
#include <algorithm>
#include <chrono>
#include <thread>

#include <Windows.h>
#include <windowsx.h>
#include <tchar.h>


#define WIN_NAME "Glass Path Tracing Window"
#define WinW_Padding 100
#define WinH_Padding 100
#define WinW_Default 400
#define WinH_Default 400
#define WinW_Min 100
#define WinH_Min 100
#define WinW_Max 800
#define WinH_Max 800

void Init();  // called before window is created
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
	case WM_GETMINMAXINFO: { LPMINMAXINFO lpMMI = (LPMINMAXINFO)lParam; lpMMI->ptMinTrackSize.x = WinW_Min, lpMMI->ptMinTrackSize.y = WinH_Min, lpMMI->ptMaxTrackSize.x = WinW_Max, lpMMI->ptMaxTrackSize.y = WinH_Max; break; }
	case WM_PAINT: { PAINTSTRUCT ps; HDC hdc = BeginPaint(hWnd, &ps), HMem = CreateCompatibleDC(hdc); HBITMAP hbmOld = (HBITMAP)SelectObject(HMem, _HIMG); BitBlt(hdc, 0, 0, _WIN_W, _WIN_H, HMem, 0, 0, SRCCOPY); SelectObject(HMem, hbmOld); EndPaint(hWnd, &ps); DeleteDC(HMem), DeleteDC(hdc); break; }
#define _USER_FUNC_PARAMS GET_X_LPARAM(lParam), _WIN_H - 1 - GET_Y_LPARAM(lParam)
	case WM_MOUSEMOVE: { MouseMove(_USER_FUNC_PARAMS); _RDBK } case WM_MOUSEWHEEL: { MouseWheel(GET_WHEEL_DELTA_WPARAM(wParam)); _RDBK }
	case WM_LBUTTONDOWN: { SetCapture(hWnd); MouseDownL(_USER_FUNC_PARAMS); _RDBK } case WM_LBUTTONUP: { ReleaseCapture(); MouseUpL(_USER_FUNC_PARAMS); _RDBK }
	case WM_RBUTTONDOWN: { MouseDownR(_USER_FUNC_PARAMS); _RDBK } case WM_RBUTTONUP: { MouseUpR(_USER_FUNC_PARAMS); _RDBK }
	case WM_SYSKEYDOWN:; case WM_KEYDOWN: { if (wParam >= 0x08) KeyDown(wParam); _RDBK } case WM_SYSKEYUP:; case WM_KEYUP: { if (wParam >= 0x08) KeyUp(wParam); _RDBK }
	} return DefWindowProc(hWnd, message, wParam, lParam);
}
#if 0
int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, PWSTR pCmdLine, int nCmdShow) {
#else
int main(int argc, char* argv[]) {
	Init();
	HINSTANCE hInstance = NULL; int nCmdShow = SW_RESTORE;
#endif
	WNDCLASSEX wc; wc.cbSize = sizeof(WNDCLASSEX), wc.style = 0, wc.lpfnWndProc = WndProc, wc.cbClsExtra = wc.cbWndExtra = 0, wc.hInstance = hInstance; wc.hIcon = wc.hIconSm = 0, wc.hCursor = LoadCursor(NULL, IDC_ARROW), wc.hbrBackground = CreateSolidBrush(RGB(0, 0, 0)), wc.lpszMenuName = NULL, wc.lpszClassName = _T(WIN_NAME); if (!RegisterClassEx(&wc)) return -1;
	_HWND = CreateWindow(_T(WIN_NAME), _T(WIN_NAME), WS_OVERLAPPEDWINDOW, WinW_Padding, WinH_Padding, WinW_Default, WinH_Default, NULL, NULL, hInstance, NULL); ShowWindow(_HWND, nCmdShow); UpdateWindow(_HWND);
	MSG message; while (GetMessage(&message, 0, 0, 0)) { TranslateMessage(&message); DispatchMessage(&message); } return (int)message.wParam;
}



#include "glslmath.h"
#undef texelFetch

#include "../octatree/bvh.h"
BVH *bvh_glass;
BVH *bvh_content;


namespace GLSL {

	// function overloading
	bool intersectGlass(vec3 ro, vec3 rd, float &t, float t1, vec3 &n) {
		vec3 color;
		if (intersectBVH(bvh_glass, ro, rd, t1, n, color)) {
			t = t1; return true;
		}
		else return false;
	}
	bool intersectContent(vec3 ro, vec3 rd, float &t, float t1, vec3 &n, vec3 &col) {
		if (intersectBVH(bvh_content, ro, rd, t1, n, col)) {
			t = t1; return true;
		}
		else return false;
	}

	// uniforms
	float iTime = 0.0;
	int iFrame = -1;
	vec3 iResolution = vec3(0.0, 0.0, 1.0);
	vec4 iMouse = vec4(0, 0, 0, 0);

	vec4 gl_FragCoord;

	// texture
	vec4 *iChannel0 = nullptr;
	vec4 texelFetch(const vec4 channel[], ivec2 coord, int plane) {
		int x = clamp(coord.x, 0, _WIN_W - 1);
		int y = clamp(coord.y, 0, _WIN_H - 1);
		return channel[y*_WIN_W + x];
	}

	// source code
#define __CPLUSPLUS
#include "../.glsl.cpp"

	// render to buffer
	vec4 *FrameBuffer = nullptr;
}


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
		if (GLSL::iChannel0) delete GLSL::iChannel0;
		GLSL::iChannel0 = new vec4[_WIN_W*_WIN_H];
		if (GLSL::FrameBuffer) delete GLSL::FrameBuffer;
		GLSL::FrameBuffer = new vec4[_WIN_W*_WIN_H];
	}

	// set uniform variables
	GLSL::iFrame += 1;
	GLSL::iResolution.x = (float)_WIN_W;
	GLSL::iResolution.y = (float)_WIN_H;
	GLSL::iTime = std::chrono::duration<float>(std::chrono::high_resolution_clock::now() - _TimeStart).count();
	for (int i = 0; i < _WIN_W*_WIN_H; i++) {
		GLSL::iChannel0[i] = GLSL::FrameBuffer[i];
	}

	// render frame
	auto time_0 = std::chrono::high_resolution_clock::now();
	batchRender([](int beg, int end, int step, bool* sig) {
		int WIN_SIZE = _WIN_W * _WIN_H;
		for (int k = beg; k < end; k += step) {
			int j = k / _WIN_W, i = k % _WIN_W;
			GLSL::gl_FragCoord = vec4(vec2(i, j) + 0.5f, 0.0f, 0.0f);
			GLSL::mainImage(GLSL::FrameBuffer[k], GLSL::gl_FragCoord.xy());
		}
		if (sig) *sig = true;
	}, _WIN_W*_WIN_H);
	float time_elapsed = std::chrono::duration<float>(std::chrono::high_resolution_clock::now() - time_0).count();

	// display image
	for (int j = 0; j < _WIN_H; j++) {
		for (int i = 0; i < _WIN_W; i++) {
			vec3 rgb = saturate(GLSL::FrameBuffer[j*_WIN_W + i].xyz());
			COLORREF c = (COLORREF(255.0*rgb.x) << 16) | (COLORREF(255.0*rgb.y) << 8) | (COLORREF(255.0*rgb.z));
			_WINIMG[j*_WIN_W + i] = c;
		}
	}

	// display information
	char text[1024];
	sprintf(text, "%.1fms (%.1f fps) [%dx%d]", 1000.0f * time_elapsed, 1.0f / time_elapsed, _WIN_W, _WIN_H);
	SetWindowTextA(_HWND, text);
	printf("iMouse=(%d,%d) iResolution=(%d,%d)\n", (int)GLSL::iMouse.x, (int)GLSL::iMouse.y, (int)GLSL::iResolution.x, (int)GLSL::iResolution.y);
}


// ============================================== User ==============================================

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "libraries/stb_image_write.h"

void render_save() {
	_WIN_W = 1200, _WIN_H = 1000;
	_WINIMG = new COLORREF[_WIN_W*_WIN_H];
	for (int i = 0; i < 1000; i++) {
		printf("%d\n", i);
		render();
	}
	for (int i = 0; i < _WIN_W*_WIN_H; i++) {
		COLORREF c = _WINIMG[i];
		c = ((c & 0x000000ff) << 16) | (c & 0x0000ff00) | ((c & 0x00ff0000) >> 16) | 0xff000000;
		_WINIMG[i] = c;
	}
	for (int i = 0; i < _WIN_W; i++) for (int j = 0; j < _WIN_H / 2; j++)
		std::swap(_WINIMG[j*_WIN_W + i], _WINIMG[(_WIN_H - j - 1)*_WIN_W + i]);
	stbi_write_png("D:\\glass.png", _WIN_W, _WIN_H, 4, _WINIMG, 4 * _WIN_W);
	exit(0);
}

void Init() {
	bvh_glass = loadModel("D://.ply");
	bvh_content = loadModel("D:/Homework/AVI4M/AVI4M-ISP/models/group_01_sdf.ply");
	render_save();
	_TimeStart = std::chrono::high_resolution_clock::now();
}

void WindowResize(int _oldW, int _oldH, int _W, int _H) {

	_RenderNeeded = true;
}
void WindowClose() {
}

void MouseWheel(int _DELTA) {
	_RenderNeeded = true;
}
void MouseDownL(int _X, int _Y) {
	GLSL::iMouse.zw() = vec2(_X, _Y);
	_RenderNeeded = true;
}
void MouseMove(int _X, int _Y) {
	if (GLSL::iMouse.z > 0.0f) {
		GLSL::iMouse.w = -abs(GLSL::iMouse.w);
		GLSL::iMouse.xy() = vec2(_X, _Y);
	}
	_RenderNeeded = true;
}
void MouseUpL(int _X, int _Y) {
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


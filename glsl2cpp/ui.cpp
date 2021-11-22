//**/

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
	case WM_CREATE: { if (!_HWND) Init(); break; } case WM_CLOSE: { WindowClose(); DestroyWindow(hWnd); return 0; } case WM_DESTROY: { PostQuitMessage(0); return 0; }
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
#if 0
int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, PWSTR pCmdLine, int nCmdShow) {
#else
int main(int argc, char* argv[]) {
	HINSTANCE hInstance = NULL; int nCmdShow = SW_RESTORE;
#endif
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



namespace GLSL {

#define GLSLMATH_USE_STD_MATH
#include "glslmath.h"

	float iTime = 0.0;
	int iFrame = -1;
	vec3 iResolution = vec3(0.0, 0.0, 1.0);
	vec4 iMouse = vec4(0, 0, 0, 0);

	vec4 gl_FragCoord;

#include "../.glsl.cpp"

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
		if (GLSL::FrameBuffer) delete GLSL::FrameBuffer;
		GLSL::FrameBuffer = new GLSL::vec4[_WIN_W*_WIN_H];
	}

	// set uniform variables
	GLSL::iFrame += 1;
	GLSL::iResolution.x = (float)_WIN_W;
	GLSL::iResolution.y = (float)_WIN_H;
	GLSL::iTime = std::chrono::duration<float>(std::chrono::high_resolution_clock::now() - _TimeStart).count();

	// render frame
	auto time_0 = std::chrono::high_resolution_clock::now();
	batchRender([](int beg, int end, int step, bool* sig) {
		int WIN_SIZE = _WIN_W * _WIN_H;
		for (int k = beg; k < end; k += step) {
			int j = k / _WIN_W, i = k % _WIN_W;
			GLSL::gl_FragCoord = GLSL::vec4(GLSL::vec2(i, j) + 0.5f, 0.0f, 0.0f);
			GLSL::mainImage(GLSL::FrameBuffer[k], GLSL::gl_FragCoord.xy());
		}
		if (sig) *sig = true;
	}, _WIN_W*_WIN_H);
	float time_elapsed = std::chrono::duration<float>(std::chrono::high_resolution_clock::now() - time_0).count();

	// display image
	for (int j = 0; j < _WIN_H; j++) {
		for (int i = 0; i < _WIN_W; i++) {
			GLSL::vec3 rgb = GLSL::saturate(GLSL::FrameBuffer[j*_WIN_W + i].xyz());
			COLORREF c = (COLORREF(255.0*rgb.x) << 16) | (COLORREF(255.0*rgb.y) << 8) | (COLORREF(255.0*rgb.z));
			_WINIMG[j*_WIN_W + i] = c;
		}
	}

	// the actual fps is less because of display time
	char text[1024];
	sprintf(text, "%.1fms (%.1f fps) [%dx%d]", 1000.0f * time_elapsed, 1.0f / time_elapsed, _WIN_W, _WIN_H);
	SetWindowTextA(_HWND, text);
}


// ============================================== User ==============================================

void Init() {
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
	GLSL::iMouse.zw() = GLSL::vec2(_X, _Y);
	_RenderNeeded = true;
}
void MouseMove(int _X, int _Y) {
	if (GLSL::iMouse.z > 0.0f) {
		GLSL::iMouse.w = -abs(GLSL::iMouse.w);
		GLSL::iMouse.xy() = GLSL::vec2(_X, _Y);
	}
	_RenderNeeded = true;
}
void MouseUpL(int _X, int _Y) {
	GLSL::iMouse.zw() = -GLSL::abs(GLSL::iMouse.zw());
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


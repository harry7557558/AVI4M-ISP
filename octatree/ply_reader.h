#pragma once

#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>



struct ply_triangle {
	int v[3];
	int& operator[] (int d) {
		return v[d];
	}
};



bool readPLY(FILE* fp, vec3* &Vs, ply_triangle* &Fs, int &VN, int &FN, vec3* &v_col) {
	if (fp == 0) return false;

	const int MAX_SIZE = 0x10000;
	int buffer_size = MAX_SIZE;
	int index = buffer_size;
	uint8_t buf[MAX_SIZE + 1]; buf[MAX_SIZE] = 0;
	auto next_byte = [&]()->int {
		if (index >= buffer_size) {
			if (buffer_size == MAX_SIZE)
				buffer_size = (int)fread(buf, 1, MAX_SIZE, fp), index = 0;
			else return EOF;
		}
		return (int)buf[index++];
	};
	auto ignore_whitespace = [&]()->bool {
		for (int i = 0;; i++) {
			if (index >= buffer_size) {
				if (buffer_size == MAX_SIZE)
					buffer_size = (int)fread(buf, 1, MAX_SIZE, fp), index = 0;
				else return true;
			}
			char c = buf[index++];
			if (c != ' ' && c != '\r' && c != '\n') {
				index--; return i != 0;
			}
		}
	};
	auto check_str = [&](const char* s)->bool {
		while (*s) {
			if (*(s++) != next_byte()) return false;
		}
		return true;
	};

	auto read_string = [&](char end_char)->std::string {
		std::string s;
		for (int i = 0; ; i++) {
			char c = next_byte();
			if (c == end_char) {
				index--;
				if (end_char == '\n' && s.back() == '\r') s.pop_back();
				return s;
			}
			if (c == 0 || c == EOF) throw(c);
			s.push_back(c);
		}
	};

	if (!check_str("ply")) return false;
	if (!ignore_whitespace()) return false;

	std::vector<std::string> header_lines;
	try {
		while (1) {
			std::string s = read_string('\n');
			if (s == "end_header") break;
			header_lines.push_back(s);
			ignore_whitespace();
		}
	} catch (...) {
		return false;
	}



	const int ASCII = 0, BINARY_BIG_ENDIAN = 1, BINARY_LITTLE_ENDIAN = 2;
	int format = -1;

	VN = FN = -1;

	int vl_xi = -1, vl_yi = -1, vl_zi = -1;  // xyz in vertex list
	int vl_ri = -1, vl_gi = -1, vl_bi = -1;  // rgb in vertex list
	int fl_li = -1;  // vertex_indices in face list
	std::vector<int> vertex_size_list, face_size_list;  // in bytes
	int vertex_property_index = 0, face_property_index = 0;
	std::string element_name = "";

	auto next_word = [](std::string &nm, std::string &s) {
		int first_space = (int)s.find(' ');
		if (first_space <= 0) return false;
		nm = s.substr(0, first_space);
		s = s.substr(first_space + 1, s.size() - first_space - 1);
		return true;
	};
	auto type_bytes = [](const std::string &type) {
		if (type == "uchar" || type == "char" || type == "int8" || type == "uint8") return 1;
		else if (type == "short" || type == "ushort" || type == "int16" || type == "uint16") return 2;
		else if (type == "float" || type == "float32" || type == "int" || type == "uint" || type == "int32" || type == "uint32") return 4;
		else if (type == "double" || type == "float64" || type == "int64" || type == "uint64") return 8;
		else return 0;
	};

	for (std::string s : header_lines) {
		std::string nm;
		if (!next_word(nm, s)) return false;

		if (nm == "format") {
			if (FN != -1 || VN != -1) return false;
			if (s == "ascii 1.0") format = ASCII;
			if (s == "binary_big_endian 1.0") format = BINARY_BIG_ENDIAN;
			if (s == "binary_little_endian 1.0") format = BINARY_LITTLE_ENDIAN;
		}

		else if (nm == "element") {
			if (!next_word(element_name, s)) return false;
			int d = std::stoi(s);
			if (element_name == "vertex") VN = d;
			else if (element_name == "face") {
				FN = d;
				if (VN == -1) return false;
			}
			else {
				if (FN == -1 || VN == -1) return false;
			}
		}

		else if (nm == "property") {
			if (element_name != "vertex" && element_name != "face") continue;
			std::string type;
			if (!next_word(type, s)) return false;
			if (element_name == "vertex") {
				if (type_bytes(type)) vertex_size_list.push_back(type_bytes(type));
				else return false;

				if (type == "float" || type == "float32") {
					if (s == "x") vl_xi = vertex_property_index;
					if (s == "y") vl_yi = vertex_property_index;
					if (s == "z") vl_zi = vertex_property_index;
				}
				else if (s == "x" || s == "y" || s == "z") return false;

				if (s == "red" || s == "green" || s == "blue") {
					if (s == "red") vl_ri = vertex_property_index;
					if (s == "green") vl_gi = vertex_property_index;
					if (s == "blue") vl_bi = vertex_property_index;
					if (type != "uchar" && type != "uint8") return false;
				}

				vertex_property_index += 1;
			}
			else if (element_name == "face") {
				if (type == "list") {
					if (!next_word(type, s)) return false;
					if (type != "uchar" && type != "char" && type != "int8" && type != "uint8") return false;
					face_size_list.push_back(type_bytes(type));
					if (!next_word(type, s)) return false;
					if (type != "uint" && type != "int" && type != "int32" && type != "uint32") return false;
					face_size_list.push_back(type_bytes(type));
					if (s != "vertex_indices" && s != "vertex_index") return false;
					if (FN == -1) return false;
					fl_li = face_property_index++;
				}
				else {
					if (type_bytes(type)) face_size_list.push_back(type_bytes(type));
					else return false;
					face_property_index += 1;
				}
			}
		}

		else if (nm != "comment") return false;
	}

	if (format == -1 || VN == -1 || FN == -1 || vl_xi == -1 || vl_yi == -1 || vl_zi == -1) return false;

	bool has_color = vl_ri != -1 && vl_gi != -1 && vl_bi != -1;
	if (has_color) v_col = new vec3[VN];


	Vs = new vec3[VN]; Fs = new ply_triangle[FN];

	if (format == ASCII) {
		return false;
	}

	else if (format == BINARY_LITTLE_ENDIAN) {

		auto read32 = [&]()->uint32_t {
			uint32_t x = 0;
			for (int i = 0; i < 4; i++) {
				uint8_t c = next_byte();
				x = x | (uint32_t(c) << (8 * i));
			}
			return x;
		};

		if (next_byte() != '\n') return false;

		float *fs = new float[vertex_property_index];
		uint32_t cr = 0, cg = 0, cb = 0;
		for (int i = 0; i < VN; i++) {
			for (int u = 0; u < vertex_property_index; u++) {
				if (u == vl_xi || u == vl_yi || u == vl_zi) *(uint32_t*)&fs[u] = read32();
				else if (u == vl_ri) cr = next_byte();
				else if (u == vl_gi) cg = next_byte();
				else if (u == vl_bi) cb = next_byte();
				else for (int _ = 0; _ < vertex_size_list[u]; _++) next_byte();
			}
			Vs[i].x = fs[vl_xi], Vs[i].y = fs[vl_yi], Vs[i].z = fs[vl_zi];
			if (has_color) v_col[i] = vec3(cr, cg, cb) / 255.0;
		}
		delete fs;

		for (int i = 0; i < FN; i++) {
			for (int u = 0; u < face_property_index; u++) {
				if (u != fl_li) {
					for (int _ = 0; _ < face_size_list[u]; _++) next_byte();
					continue;
				}
				int n = (int)next_byte();
				if (n != 3) return false;
				for (int u = 0; u < 3; u++) {
					int d = (int)read32();
					if (d >= 0 && d < VN) Fs[i][u] = d;
					else return false;
				}
			}
		}

	}

	else if (format == BINARY_BIG_ENDIAN) {
		return false;
	}

	else return false;


	return true;
}

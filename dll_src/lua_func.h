#pragma once

#include <lua.h>

#include "structure.h"

void getpixeldata(lua_State *L, Pixel_RGBA **pixels, int *w, int *h);
void putpixeldata(lua_State *L, Pixel_RGBA *pixels);
int error(lua_State *L, char *s);
void print_stack(lua_State *L);

#include <stdio.h>
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include "structure.h"
#include "lua_func.h"


//obj.getpixeldataをつかってピクセルを入手する関数
void getpixeldata(lua_State *L, Pixel_RGBA **pixels, int *w, int *h) {
    lua_getglobal(L, "obj");
    lua_getfield(L, -1, "getpixeldata");
    lua_call(L, 0, 3);
    *h = lua_tointeger(L, -1);
    lua_pop(L, 1);
    *w = lua_tointeger(L, -1);
    lua_pop(L, 1);
    *pixels = (Pixel_RGBA*)lua_touserdata(L, -1);
    lua_pop(L, 1);
}

//obj.putpixeldataをつかって描画する関数
void putpixeldata(lua_State *L, Pixel_RGBA *pixels) {
    lua_getfield(L, -1, "putpixeldata");
    lua_pushlightuserdata(L, pixels);
    lua_call(L, 1, 0);
}

//エラー処理
int error(lua_State *L, char *s) {
    fprintf(stderr, "\033[31m%s\033[0m\n", s);
    lua_pushboolean(L, 0);
    lua_pushstring(L, s);
    return 2;
}

//luaのスタックの様子を表示する関数(デバッグ用)
void print_stack(lua_State *L) {
    int top = lua_gettop(L);
    printf("スタックの内容（合計 %d 個）:\n", top);
    for (int i = 1; i <= top; i++) {
        int t = lua_type(L, i);
        printf("Stack[%2d-%10s] : ", i, lua_typename(L, t));
        switch (t) {
        case LUA_TSTRING:
            printf("%d: '%s'\n", i, lua_tostring(L, i));
            break;
        case LUA_TBOOLEAN:
            printf("%d: %s\n", i, lua_toboolean(L, i) ? "true" : "false");
            break;
        case LUA_TNUMBER:
            printf("%d: %g\n", i, lua_tonumber(L, i));
            break;
        default:
            printf("%d: %s\n", i, lua_typename(L, t));
            break;
        }
    }
}

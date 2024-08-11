//ISO C17
#include <stdio.h>
#include <stdlib.h>
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
#include <omp.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "structure.h"
#include "lua_func.h"

#define MIN(i, j) (((i) < (j)) ? i : j)
#define MAX(i, j) (((i) > (j)) ? i : j)
#define CLAMP(x, min, max) (((x) <= min) ? min : (((x) >= max) ? max : x))


static void angle_to_rg(double angle, double intensity, unsigned char *r, unsigned char *g);
static int check_color(int len, int *list, int mode, Pixel_RGBA *pixels);
static double calculate_target(int w, int h, int subdivison_target, double subdivison_target_scale);
static inline int is_inside_regular_polygon(double norm_dot, double radius, double initial_radius, int map_pattern, int normal_divisions);
static inline void calculate_intensity(double dx, double dy, double target, int is_normal, Vec2 *std_vec, Param *param, double *intensity, int *n);
static int draw_polygon(lua_State *L, Pixel_RGBA *pixels, int w, int h, const Param *param);
static int draw_rectangle(lua_State *L, Pixel_RGBA *pixels, int w, int h, Param *param);
//未実装
//static int draw_rotation(Pixel_RGBA *pixels, int w, int h, Param *param);


int edit(lua_State *L) {
    //パラメータを構造体で管理
    //ここで条件外の数値を弾く
    Param param = {0};
    param.mode = lua_isnumber(L, 1) ? lua_tointeger(L, 1) : 1;
    if (param.mode < 1 || param.mode > 3)
        return error(L, "[DisplacementMapEditor.dll] This mode doesn't exist.");

    param.direction = lua_isnumber(L, 2) ? M_PI * lua_tonumber(L, 2) / 180.0 : 0.0; // deg to rad
    param.intensity = lua_isnumber(L, 3) ? lua_tonumber(L, 3) : 1.0;
    if (param.intensity < 0.0)
        return error(L, "[DisplacementMapEditor.dll] This intensity is out of range.");

    param.initial_angle = lua_isnumber(L, 4) ? M_PI * lua_tonumber(L, 4) / 180.0 : 0.0; // deg to rad
    param.divisions = lua_isnumber(L, 5) ? lua_tointeger(L, 5) : 0;
    if (param.divisions < 1)
        return error(L, "[DisplacementMapEditor.dll] This number of divisions is out of range.");

    param.color_selection_mode = lua_isnumber(L, 6) ? lua_tointeger(L, 6) : 1;
    if (param.color_selection_mode < 1 || param.color_selection_mode > 3)
        return error(L, "[DisplacementMapEditor.dll] This color selection mode doesn't exist.");

    if (lua_istable(L, 7))
        param.color_table_index = 7;
    else
        return error(L, "[DisplacementMapEditor.dll] This parameter is not table.");

    param.flip = lua_isboolean(L, 8) ? lua_toboolean(L, 8) : 0;
    param.map_type = lua_isnumber(L, 9) ?  lua_tointeger(L, 9) : 1;
    if (param.map_type < 1 || param.map_type > 2)
        return error(L, "[DisplacementMapEditor.dll] This map type doesn't exist.");

    param.map_pattern = lua_isnumber(L, 10) ? lua_tointeger(L, 10) : 1;
    if (param.map_pattern < 1 || param.map_pattern > 3)
        return error(L, "[DisplacementMapEditor.dll] This map pattern doesn't exist.");

    if (lua_istable(L, 11))
        param.subdivison_target_table_index = 11;
    else
        return error(L, "[DisplacementMapEditor.dll] This parameter is not table.");
    
    if (lua_objlen(L, param.subdivison_target_table_index) < 1 || lua_objlen(L, param.subdivison_target_table_index) > 2)
        return error(L, "[DisplacementMapEditor.dll] The length of this subdivison target table is out of range.");

    if (lua_istable(L, 12))
        param.subdivison_target_scale_table_index = 12;
    else
        return error(L, "[DisplacementMapEditor.dll] This parameter is not table.");

    if (lua_objlen(L, param.subdivison_target_scale_table_index) < 1 || lua_objlen(L, param.subdivison_target_scale_table_index) > 2)
        return error(L, "[DisplacementMapEditor.dll] The length of this subdivison target scale table is out of range.");

    param.normal_divisions = lua_isnumber(L, 13) ? lua_tointeger(L, 13) : 0;
    if (param.normal_divisions < 0)
        return error(L, "[DisplacementMapEditor.dll] This normal divisions is out of range.");

    int max_threads = omp_get_max_threads();
    param.threads = lua_isnumber(L, 14) ? CLAMP(lua_tointeger(L, 14), 0, max_threads) : max_threads;
    if (param.threads == 0)
        param.threads = max_threads;

    //ピクセルデータを入手
    Pixel_RGBA *pixels;
    int w, h;

    getpixeldata(L, &pixels, &w, &h);

    
    switch (param.mode) {
        //彩色mode
        case 1: {
            unsigned char r, g;
            angle_to_rg(param.direction, param.intensity, &r, &g);

            int i = 0;
            #pragma omp parallel for num_threads(param.threads)
            for (i = 0; i < w * h; i++) {
                pixels[i].r = r;
                pixels[i].g = g;
                pixels[i].b = 128;
            }

            break;
        }
        //修正mode
        case 2: {
            int color_table_len = lua_objlen(L, param.color_table_index);
            int *color_list = (int*)calloc(color_table_len, sizeof(int));
            if (!color_list) {
                color_table_len = 0;
                free(color_list);
                return error(L, "[DisplacementMapEditor.dll] Memory allocation failed.");
            }
            for (int i = 1; i <= color_table_len; i++) {
                lua_pushinteger(L, i);
                lua_gettable(L, param.color_table_index);
                if (lua_isnumber(L, -1))
                    color_list[i - 1] = lua_tointeger(L, -1);
                else
                    fprintf(stderr, "[DisplacementMapEditor.dll] %d番目の要素が数字ではありません。\n", i);

                lua_pop(L, 1);
            }

            double UVEC_X = cos(param.direction), UVEC_Y = sin(param.direction); //the x and y components of unit vector

            int i = 0;
            #pragma omp parallel for num_threads(param.threads)
            for (i = 0; i < w * h; i++) {
                if ((pixels[i].a == 0) || check_color(color_table_len, color_list, param.color_selection_mode, (pixels + i)))
                    continue;

                unsigned char R = pixels[i].r, G = pixels[i].g;
                pixels[i].r = (unsigned char)CLAMP(ceil(param.intensity * ((R - 128) * UVEC_X - (G - 128) * UVEC_Y) + 128), 0, 255);
                pixels[i].g = (unsigned char)CLAMP(ceil(param.intensity * ((G - 128) * UVEC_X + (R - 128) * UVEC_Y) + 128), 0, 255);
            }

            free(color_list);

            break;
        }
        //map生成mode
        case 3: {
            int check = 0;
            switch (param.map_type) {
                case 1:
                    check = draw_polygon(L, pixels, w, h, &param);
                    break;
                case 2:
                    check = draw_rectangle(L, pixels, w, h, &param);
                    break;
                default:
                    break;
            }
            if (check == -1)
                return 2;

            break;
        }
        default:
            break;
    }

    //描画
    putpixeldata(L, pixels);
    lua_pushboolean(L, 1);
    lua_pushstring(L, "[DisplacementMapEditor.dll] success");

    return 2;
}

//角度から色(RG)を求める関数
static void angle_to_rg(double angle, double intensity, unsigned char *r, unsigned char *g) {
    *r = (unsigned char)ceil(127.5 * (1 - MIN(intensity, 1) * cos(angle)));
    *g = (unsigned char)ceil(127.5 * (1 - MIN(intensity, 1) * sin(angle)));
}

//テーブルに色が含まれているかどうかチェックする関数
static int check_color(int len, int *list, int mode, Pixel_RGBA *pixels) {
    switch (mode) {
        case 1:
            return 0;
            break;
        case 2:
        case 3:
            int hexcolor = (pixels->r << 16) | (pixels->g << 8) | pixels->b;
            for (int j = 0; j < len; j++)
                if (hexcolor == list[j])
                    return (mode == 2) ? 1 : 0;

            return (mode == 2) ? 0 : 1;
            break;
        default:
            return 1;
            break;
    }
}

//分割ターゲットを求める関数
static double calculate_target(int w, int h, int subdivison_target, double subdivison_target_scale) {
    double target = 0.0;
    switch (subdivison_target) {
    case 1:
        target = subdivison_target_scale * MIN(w, h);
        break;
    case 2:
        target = subdivison_target_scale * MAX(w, h);
        break;
    case 3:
        target = subdivison_target_scale * sqrt(w * w + h * h);
        break;
    case 4:
        target = subdivison_target_scale;
        break;
    default:
        return 0.0;
    }

    return target;
}

//正n角形内部に頂点が含まれているかチェックする関数
static inline int is_inside_regular_polygon(double norm_dot, double radius, double initial_radius, int map_pattern, int normal_divisions) {
    if (radius == 0 || normal_divisions == 0)
        return 1;

    if (norm_dot <= radius)
        return (normal_divisions + 1) % 2;

    if (normal_divisions > 1) {
        switch (map_pattern) {
        case 1:
            radius += initial_radius;
            break;
        case 2:
            radius *= 2.0;
            break;
        }
        return is_inside_regular_polygon(norm_dot, radius, initial_radius, map_pattern, normal_divisions - 1);
    }

    return normal_divisions % 2;
}

//座標から強度を計算する関数
static inline void calculate_intensity(double dx, double dy, double target, int is_normal, Vec2 *std_vec, Param *param, double *intensity, int *section) {
    int divisions = is_normal ? param->normal_divisions : param->divisions;
    double width = target / (divisions + 1.0);

    if (divisions == 0 || width == 0)
        *intensity = 2.0 * MIN(param->intensity, 1.0);

    double d = is_normal ? fabs(std_vec->y * dx - std_vec->x * dy) : fabs(std_vec->x * dx + std_vec->y * dy);
    double dot = dx * std_vec->x + dy * std_vec->y;
    double length = (divisions % 2) ? d : d - 0.5 * width;
    double pre_intensity = 0.0;
    if (!(divisions % 2) && length < 0.0) {
        pre_intensity = 0.0;
    } else if (divisions == 1 || divisions == 2) {
        pre_intensity = 1.0;
    } else {
        int n = MIN((int)(length / width), (divisions - 1) / 2);
        int range = divisions % 2 ? (divisions - 1) / 2 : (divisions - 2) / 2;
        pre_intensity = n * (2.0 * MIN(param->intensity, 1.0) - 1.0) / (double)range + 1.0;
        *section = n;
    }

    if (dot >= 0)
        *intensity = pre_intensity;
    else
        *intensity = -pre_intensity;
}

//正n角形状のmapを生成する関数
static int draw_polygon(lua_State *L, Pixel_RGBA *pixels, int w, int h, const Param *param) {
    double cx = w / 2.0, cy = h / 2.0;
    double half_fan_angle = M_PI / param->divisions;
    double range = cos(half_fan_angle);
    lua_pushinteger(L, 1);
    lua_gettable(L, param->subdivison_target_table_index);
    int subdivison_target = lua_isnumber(L, -1) ? lua_tointeger(L, -1) : 0;
    lua_pop(L, 1);
    if (subdivison_target < 1 || subdivison_target > 4) {
        error(L, "[DisplacementMapEditor.dll] This division target doesn't exist.");
        return -1;
    }

    lua_pushinteger(L, 1);
    lua_gettable(L, param->subdivison_target_scale_table_index);
    double subdivison_target_scale = lua_isnumber(L, -1) ? lua_tonumber(L, -1) : -1.0;
    lua_pop(L, 1);
    if (subdivison_target < 0.0) {
        error(L, "[DisplacementMapEditor.dll] This scale of this division target is out of range.");
        return -1;
    }

    double radius = calculate_target(w, h, subdivison_target, subdivison_target_scale);
    double radius_range = cos(M_PI / param->divisions);
    switch (param->map_pattern) {
        case 1:
            radius = radius * radius_range / (param->normal_divisions + 1);
            break;
        case 2:
            radius = radius * radius_range / pow(2.0, param->normal_divisions);
            break;
        default:
            return -1;
    }
    Vec2 rot = { cos(2.0 * half_fan_angle), sin(2.0 * half_fan_angle) };
    Vec2 initial_vec = { cos(half_fan_angle + param->initial_angle), sin(half_fan_angle + param->initial_angle) };
    Pixel_RG *pixel_rg = (Pixel_RG*)calloc(param->divisions, sizeof(Pixel_RG));
    if (!pixel_rg) {
        free(pixel_rg);
        error(L, "[DisplacementMapEditor.dll] Memory allocation failed.");
        return -1;
    }

    for (int i = 0; i < param->divisions; i++) {
        unsigned char r, g;
        angle_to_rg(half_fan_angle * (1 + 2.0 * i) + param->initial_angle + param->direction, param->intensity, &r, &g);
        pixel_rg[i].r = r;
        pixel_rg[i].g = g;
    }

    int i = 0;
    #pragma omp parallel for num_threads(param->threads)
    for (i = 0; i < w * h; i++) {
        int x = i % w, y = i / w;
        double dx = x - cx, dy = y - cy;
        double norm = sqrt(dx * dx + dy * dy);
        if (norm == 0)
            norm = 1e-10;

        Vec2 std_vec = initial_vec;
        double dot = 0;
        int section = 0;
        for ( ; section < param->divisions; section++) {
            dot = (dx * std_vec.x + dy * std_vec.y) / norm;
            if (dot >= range)
                break;

            std_vec = (Vec2){
                std_vec.x * rot.x - std_vec.y * rot.y,
                std_vec.x * rot.y + std_vec.y * rot.x
            };
        }

        int col_flip = (is_inside_regular_polygon(norm * dot, radius, radius, param->map_pattern, param->normal_divisions) != (param->flip && (section % 2))) ? 1 : 0;

        pixels[i].r = col_flip ? pixel_rg[section].r : 255 - pixel_rg[section].r;
        pixels[i].g = col_flip ? pixel_rg[section].g : 255 - pixel_rg[section].g;
        pixels[i].b = 128;
    }

    free(pixel_rg);
    return 0;
}

//長方形状のマップを生成する関数
static int draw_rectangle(lua_State *L, Pixel_RGBA *pixels, int w, int h, Param *param) {
    double cx = w / 2.0, cy = h / 2.0;
    int subdivison_target_table[2] = { 0 };
    int subdivison_target_table_len = lua_objlen(L, param->subdivison_target_table_index);
    for (int i = 1; i <= subdivison_target_table_len; i++) {
        lua_pushinteger(L, i);
        lua_gettable(L, param->subdivison_target_table_index);
        subdivison_target_table[i - 1] = lua_isnumber(L, -1) ? lua_tointeger(L, -1) : 0;
        lua_pop(L, 1);
        if (subdivison_target_table[i - 1] < 1 || subdivison_target_table[i - 1] > 4) {
            error(L, "[DisplacementMapEditor.dll] This division target doesn't exist.");
            return -1;
        }
    }

    double subdivison_target_scale_table[2] = { 0 };
    int subdivison_target_scale_table_len = lua_objlen(L, param->subdivison_target_scale_table_index);
    for (int i = 1; i <= subdivison_target_scale_table_len; i++) {
        lua_pushinteger(L, i);
        lua_gettable(L, param->subdivison_target_scale_table_index);
        subdivison_target_scale_table[i - 1] = lua_isnumber(L, -1) ? lua_tonumber(L, -1) : -1.0;
        lua_pop(L, 1);
        if (subdivison_target_scale_table[i - 1] < 0.0) {
            error(L, "[DisplacementMapEditor.dll] This scale of this division target is out of range.");
            return -1;
        }
    }

    double target_x, target_y;
    if (subdivison_target_table_len == 1) {
        if (subdivison_target_scale_table_len == 1) {
            target_x = calculate_target(w, h, subdivison_target_table[0], subdivison_target_scale_table[0]);
            target_y = target_x;
        } else {
            target_x = calculate_target(w, h, subdivison_target_table[0], subdivison_target_scale_table[0]);
            target_y = calculate_target(w, h, subdivison_target_table[0], subdivison_target_scale_table[1]);
        }
    } else {
        if (subdivison_target_scale_table_len == 1) {
            target_x = calculate_target(w, h, subdivison_target_table[0], subdivison_target_scale_table[0]);
            target_y = calculate_target(w, h, subdivison_target_table[1], subdivison_target_scale_table[0]);
        }
        else {
            target_x = calculate_target(w, h, subdivison_target_table[0], subdivison_target_scale_table[0]);
            target_y = calculate_target(w, h, subdivison_target_table[1], subdivison_target_scale_table[1]);
        }
    }

    Vec2 std_vec[2] = {0};
    std_vec[0] = (Vec2){
        cos(param->initial_angle),
        sin(param->initial_angle)
    };
    if (param->normal_divisions != 0) {
        std_vec[1] = (Vec2){
            cos(param->initial_angle + 0.75 * M_PI),
            sin(param->initial_angle + 0.75 * M_PI)
        };
    } else {
        std_vec[1] = std_vec[0];
    }
    Pixel_RG pixel_rg[2] = {0};
    unsigned char r, g;
    if (param->normal_divisions == 0) {
        angle_to_rg(param->initial_angle + param->direction, 0.5, &r, &g);
        pixel_rg[0].r = r;
        pixel_rg[0].g = g;
    } else {
        for (int i = 0; i < 2; i++) {
            angle_to_rg(param->initial_angle + param->direction + M_PI * (1.0 - 2.0 * i) / 4.0, 0.5, &r, &g);
            pixel_rg[i].r = r;
            pixel_rg[i].g = g;
        }
    }
    int i = 0;
    #pragma omp parallel for num_threads(param->threads)
    for (i = 0; i < w * h; i++) {
        int x = i % w, y = i / w;
        double dx = x - cx, dy = y - cy;
        int col_idx = 0;
        if (param->normal_divisions != 0) {
            double norm = sqrt(dx * dx + dy * dy);
            if (norm == 0)
                norm = 1e-10;

            double dot = sqrt(2) * (dx * std_vec[1].x + dy * std_vec[1].y) / norm;
            if (dot > 1.0 || dot < -1.0)
                col_idx = 1;
        }

        double intensity = 0.0;
        int section = 0;
        int section_x = 0, section_y = 0;
        if (param->normal_divisions == 0) {
            calculate_intensity(dx, dy, target_x, 0, &std_vec[0], param, &intensity, &section);
        } else {
            double intensity_x = 0.0, intensity_y = 0.0;
            calculate_intensity(dx, dy, target_x, 0, &std_vec[0], param, &intensity_x, &section_x);
            calculate_intensity(dx, dy, target_y, 1, &std_vec[0], param, &intensity_y, &section_y);
            if ((intensity_x == 0) || (intensity_y == 0))
                intensity = 0.0;
            else if (intensity_x > 0)
                intensity = MAX(intensity_x, intensity_y);
            else {
                intensity = MIN(intensity_x, intensity_y);
            }
        }

        int col_flip = 0;
        if (param->normal_divisions == 0)
            col_flip = param->flip && (section % 2) ? -1 : 1;
        else
            col_flip = param->flip && (MAX(section_x, section_y) % 2) ? -1 : 1;

        pixels[i].r = (unsigned char)CLAMP(col_flip * intensity * (pixel_rg[col_idx].r - 128) + 128, 0, 255);
        pixels[i].g = (unsigned char)CLAMP(col_flip * intensity * (pixel_rg[col_idx].g - 128) + 128, 0, 255);
        pixels[i].b = 128;
    }

    return 0;
}
/*
//未実装
static int draw_rotation(Pixel_RGBA *pixels, int w, int h, Param *param) {
    int i = 0;
    #pragma omp parallel for num_threads(param->threads)
    for (i = 0; i < w * h; i++) {
        int x = i % w, y = i / w;
        //cs求める
        double width = (double)w / (param->divisions + 1.0);
        double height = (double)h / (param->normal_divisions + 1.0);
        double dx = x - 0.5 * w;
        double dy = y - 0.5 * h;
        int nx = (int)(dx / width);
        int ny = (int)(dy / height);
        //double cx = (dx > 0 ? 0.5 : -0.5) * width + width * nx;
        //double cy = (dy > 0 ? 0.5 : -0.5) * height + height * ny;
        double cx = 0, cy = 0;

        unsigned char r, g;
        angle_to_rg(atan2(dy - cy, dx - cx) + 0.5 * M_PI + param->direction, param->intensity, &r, &g);
        pixels[i].r = r;
        pixels[i].g = g;
        pixels[i].b = 128;
    }

    return 0;
}
*/


static const struct luaL_reg functions[] = {
    { "edit", edit },
    { NULL, NULL }
};

__declspec(dllexport) int luaopen_DisplacementMapEditor(lua_State *L) {
    luaL_register(L, "DisplacementMapEditor", functions);
    return 1;
}

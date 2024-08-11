#pragma once

typedef struct {
    unsigned char b, g, r, a;
} Pixel_RGBA;

typedef struct {
    unsigned char g, r;
} Pixel_RG;

typedef struct {
    double x, y;
} Vec2;

typedef struct {
    int mode;
    double direction;
    double intensity;
    double initial_angle;
    int divisions;
    int color_selection_mode;
    int color_table_index;
    int flip;
    int map_type;
    int map_pattern;
    int subdivison_target_table_index;
    int subdivison_target_scale_table_index;
    int normal_divisions;
    int threads;
} Param;

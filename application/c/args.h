#ifndef ARGS_H_
#define ARGS_H_

#include "aether/aether.h"

typedef u8 OptKind;
enum OptKind_
{
    OptKind_Flag,
    OptKind_U16,
    OptKind_U32,
    OptKind_F64,
    OptKind_Str
};

typedef struct OptSpec
{
    str8    name;
    str8    alias;
    OptKind kind;
    void*   out;
    str8    metavar;
    str8    help;
    b8      positional;
} OptSpec;

static inline OptSpec* find_opt(OptSpec* opts, u64 count, str8 key)
{
    if (!opts) return NULL;

    for (u64 i = 0; i < count; ++i)
    {
        if (str8_eq(opts[i].name, key)) return &opts[i];
        if (str8_eq(opts[i].alias, key) && opts[i].alias.size) return &opts[i];
    }
    return NULL;
}


static void print_help(OptSpec* opts, u64 count)
{

    printf("\nUSAGE\n");
    printf("  laplace.exe");
    for (u64 i = 0; i < count; ++i)
    {
        const OptSpec* o = &opts[i];
        if (!o->positional) continue;
        printf(" <" STR8_FMT ">", STR8_ARG(o->name));
    }
    for (u64 i = 0; i < count; ++i)
    {
        const OptSpec* o = &opts[i];
         if (o->positional) continue;
        printf(" [");
        printf(STR8_FMT, STR8_ARG(o->name));
        if (o->kind != OptKind_Flag) printf(" <" STR8_FMT ">", STR8_ARG(o->metavar));
        printf("]");
    }
    printf("\n");

    printf("\nARGUMENTS\n");
    for (u64 i = 0; i < count; ++i)
    {
        const OptSpec* o = &opts[i];
        if (!o->positional) continue;
        printf(" " STR8_FMT "\t\t" STR8_FMT " \n", STR8_ARG(o->name), STR8_ARG(o->help));
    }

    printf("\nFLAGS\n");
    for (u64 i = 0; i < count; ++i)
    {
        const OptSpec* o = &opts[i];
        if (o->positional || o->kind != OptKind_Flag) continue;
        if (o->alias.size > 0) printf(" " STR8_FMT ", ", STR8_ARG(o->alias));
        else                   printf("     ");
        printf(STR8_FMT "\t\t" STR8_FMT "\n",
               STR8_ARG(o->name), STR8_ARG(o->help));
    }

    printf("\nOPTIONS\n");
    for (u64 i = 0; i < count; ++i)
    {
        const OptSpec* o = &opts[i];
        if (o->positional || o->kind == OptKind_Flag) continue;
        if (o->alias.size > 0) printf(" " STR8_FMT ", ", STR8_ARG(o->alias));
        else                   printf("     ");
        printf(STR8_FMT " <" STR8_FMT ">\t" STR8_FMT "\n",
               STR8_ARG(o->name), STR8_ARG(o->metavar), STR8_ARG(o->help));
    }
}

static inline b8 parse_args(int argc, char** argv, OptSpec* opts, u64 count)
{

    u64 next_pos = 0;

    for (int i = 1; i < argc; ++i)
    {
        str8 arg = str8_from_c_str(argv[i]);

        /* positional: match the Nth positional-marked OptSpec entry, in table order */
        if (arg.size == 0 || arg.data[0] != '-')
        {
            OptSpec* opt = NULL;
            for (u64 k = 0, seen = 0; k < count; ++k)
            {
                if (!opts[k].positional) continue;
                if (seen == next_pos) { opt = &opts[k]; break; }
                ++seen;
            }
            if (!opt)
            {
                fprintf(stderr, "unexpected positional argument: " STR8_FMT "\n", STR8_ARG(arg));
                return false;
            }
            *(str8*)opt->out = arg;
            ++next_pos;
            continue;
        }

        /* keyword: --name / --name=vale / --name vale, matched by name or alias */
        str8 key, inline_val;
        b8 has_inline = str8_cut(arg, STR("="), &key, &inline_val);

        OptSpec* opt = find_opt(opts, count, key);
        if (!opt) { fprintf(stderr, "Unrecognized option: " STR8_FMT "\n", STR8_ARG(key)); return false; }

        str8 val = {0};

        if (opt->kind != OptKind_Flag)
        {
            if (has_inline) val = inline_val;
            else if (i + 1 < argc) val = str8_from_c_str(argv[++i]);
            else { fprintf(stderr, STR8_FMT " requires a value\n", STR8_ARG(key)); return false; }
        }

        switch (opt->kind)
        {
            case OptKind_Flag: *(b8*)opt->out = true; break;
            case OptKind_U16: {
                u64 v;
                if (!str8_to_u64(val, &v) || v == 0 || v > U16_MAX) {
                    fprintf(stderr, STR8_FMT " value must fit in a u16, passed \n" STR8_FMT, STR8_ARG(key), STR8_ARG(val));
                    return false;
                }
                *(u16*)opt->out = (u16)v;
            } break;
            case OptKind_U32:{
                u64 v;
                if (!str8_to_u64(val, &v) || v == 0 || v > U32_MAX) {
                    fprintf(stderr, STR8_FMT " value must fit in a u32, passed \n" STR8_FMT, STR8_ARG(key), STR8_ARG(val));
                    return false;
                }
                *(u32*)opt->out = (u32)v;
            } break;
            case OptKind_F64: {
                f64 v; 
                if (!str8_to_f64(val, &v))
                {
                    fprintf(stderr, STR8_FMT " value must be a valid f64, passed \n" STR8_FMT, STR8_ARG(key), STR8_ARG(val));
                    return false;
                }
                *(f64*)opt->out = v;
            } break;
            case OptKind_Str: {
                *(str8*)opt->out = val;
            } break;

            default: {
                fprintf(stderr, "Unknown option kind\n");
                return false;
            }
        }
    }
    return true;
}


#endif // ARGS_H_

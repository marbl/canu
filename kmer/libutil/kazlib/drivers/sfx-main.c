#include <stdlib.h>

int main(int argc, char **argv)
{
    char expr_buf[256];
    char *expr, *ptr;
    sfx_rating_t eff;

    for (;;) {
	if (argc < 2) {
	    expr = expr_buf;
	    if (fgets(expr_buf, sizeof expr_buf, stdin) == 0)
		break;
	    if ((ptr = strchr(expr_buf, '\n')) != 0)
		*ptr = 0;
	} else {
	    expr = (argv++)[1];
	    if (!expr)
		break;
	}

	if (!sfx_determine(expr, &eff)) {
	    printf("expression '%s' has a syntax error\n", expr);
	    return EXIT_FAILURE;
	}

	switch (eff) {
	case sfx_none:
	    printf("expression '%s' has no side effects\n", expr);
	    break;
	case sfx_potential:
	    printf("expression '%s' may have side effects\n", expr);
	    break;
	case sfx_certain:
	    printf("expression '%s' has side effects\n", expr);
	    break;
	}
    }

    return 0;
}

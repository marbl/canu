#include <stdio.h>
#include <ctype.h>

static void cleanup(void *arg)
{
    printf("cleanup(\"%s\") called\n", (char *) arg);
}

static void bottom_level(void)
{
    char buf[256];
    printf("throw exception? "); fflush(stdout);
    fgets(buf, sizeof buf, stdin);

    if (buf[0] >= 0 && toupper(buf[0]) == 'Y')
	except_throw(1, 1, "nasty exception");
}

static void top_level(void)
{
    except_cleanup_push(cleanup, "argument");
    bottom_level();
    except_cleanup_pop(0);
}

int main(int argc, char **argv)
{
    static const except_id_t catch[] = { { 1, 1 }, { 1, 2 } };
    except_t *ex;

    /*
     * Nested exception ``try blocks''
     */

    /* outer */
    except_try_push(catch, 2, &ex);
    if (!ex) {
	/* inner */
	except_try_push(catch, 2, &ex);
	if (!ex) {
	    top_level();
	} else {
	    /* inner catch */
	    printf("caught exception (inner): \"%s\", s=%ld, c=%ld\n",
		    except_message(ex), except_group(ex), except_code(ex));
	    except_rethrow(ex);
	}
	except_try_pop();
    } else {
	/* outer catch */
	printf("caught exception (outer): \"%s\", s=%ld, c=%ld\n",
		except_message(ex), except_group(ex), except_code(ex));
    }
    except_try_pop();
    except_throw(99, 99, "exception in main");
    return 0;
}

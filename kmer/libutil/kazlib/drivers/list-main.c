#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>

typedef char input_t[256];

static int tokenize(char *string, ...)
{
    char **tokptr; 
    va_list arglist;
    int tokcount = 0;

    va_start(arglist, string);
    tokptr = va_arg(arglist, char **);
    while (tokptr) {
	while (*string && isspace((unsigned char) *string))
	    string++;
	if (!*string)
	    break;
	*tokptr = string;
	while (*string && !isspace((unsigned char) *string))
	    string++;
	tokptr = va_arg(arglist, char **);
	tokcount++;
	if (!*string)
	    break;
	*string++ = 0;
    }
    va_end(arglist);

    return tokcount;
}

static int comparef(const void *key1, const void *key2)
{
    return strcmp(key1, key2);
}

static char *dupstring(char *str)
{
    int sz = strlen(str) + 1;
    char *new = malloc(sz);
    if (new)
	memcpy(new, str, sz);
    return new;
}

int main(void)
{
    input_t in;
    list_t *l = list_create(LISTCOUNT_T_MAX);
    lnode_t *ln;
    char *tok1, *val;
    int prompt = 0;

    char *help =
	"a <val>                append value to list\n"
	"d <val>                delete value from list\n"
	"l <val>                lookup value in list\n"
	"s                      sort list\n"
	"c                      show number of entries\n"
	"t                      dump whole list\n"
	"p                      turn prompt on\n"
	"q                      quit";

    if (!l)
	puts("list_create failed");

    for (;;) {
	if (prompt)
	    putchar('>');
	fflush(stdout);

	if (!fgets(in, sizeof(input_t), stdin))
	    break;

	switch(in[0]) {
	    case '?':
		puts(help);
		break;
	    case 'a':
		if (tokenize(in+1, &tok1, (char **) 0) != 1) {
		    puts("what?");
		    break;
		}
		val = dupstring(tok1);
		ln = lnode_create(val);
	
		if (!val || !ln) {
		    puts("allocation failure");
		    if (ln)
			lnode_destroy(ln);
		    free(val);
		    break;
		}
    
		list_append(l, ln);
		break;
	    case 'd':
		if (tokenize(in+1, &tok1, (char **) 0) != 1) {
		    puts("what?");
		    break;
		}
		ln = list_find(l, tok1, comparef);
		if (!ln) {
		    puts("list_find failed");
		    break;
		}
		list_delete(l, ln);
		val = lnode_get(ln);
		lnode_destroy(ln);
		free(val);
		break;
	    case 'l':
		if (tokenize(in+1, &tok1, (char **) 0) != 1) {
		    puts("what?");
		    break;
		}
		ln = list_find(l, tok1, comparef);
		if (!ln)
		    puts("list_find failed");
		else
		    puts("found");
		break;
	    case 's':
		list_sort(l, comparef);
		break;
	    case 'c':
		printf("%lu\n", (unsigned long) list_count(l));
		break;
	    case 't':
		for (ln = list_first(l); ln != 0; ln = list_next(l, ln))
		    puts(lnode_get(ln));
		break;
	    case 'q':
		exit(0);
		break;
	    case '\0':
		break;
	    case 'p':
		prompt = 1;
		break;
	    default:
		putchar('?');
		putchar('\n');
		break;
	}
    }

    return 0;
}

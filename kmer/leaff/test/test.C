#include "fastafile.H"

char *filename = 0L;

FILE *openTmpFile(void) {
  FILE *O = fopen("test.C.tmp.out", "w");
  if (O == 0L) {
    fprintf(stderr, "can't open test.C.tmp.out\n");
    exit(1);
  }
  return(O);
}

void closeTmpFile(FILE *O) {
  fclose(O);
}

void removeTmpFile(void) {
  unlink("test.C.tmp.out");
}





void
test1(int style) {
  FastASequenceInCore  *x;
  FILE *O = openTmpFile();

  FastAFileWrapper  F(filename);

  switch (style) {
    case 1:
      //  nop!
      break;
    case 2:
      F.openIndex();
      break;
    case 3:
      F.setIndexToSaveIDs();
      F.openIndex();
      break;
    case 4:
      F.setIndexToSaveDeflines();
      F.openIndex();
      break;
  }

  while (!F.eof()) {
    x = F.getSequence();
    fprintf(O, "%s\n%s\n", x->header(), x->sequence());
    delete x;
  }

  closeTmpFile(O);
  removeTmpFile();
}



void
test2(int style, char *target) {
  FastASequenceInCore  *x;

  FastAFileWrapper  F(filename);

  switch (style) {
    case 5:
      F.openIndex();
      if (F.find(atoi(target)) == false) {
        fprintf(stderr, "Failed to find target = %d for test 5\n", atoi(target));
        exit(1);
      }
      break;
    case 6:
      F.setIndexToSaveIDs();
      F.openIndex();
      if (F.find(target) == false) {
        fprintf(stderr, "Failed to find target = '%s' for test 6\n", target);
        exit(1);
      }
      break;
    case 7:
      F.setIndexToSaveDeflines();
      F.openIndex();
      if (F.find(target) == false) {
        fprintf(stderr, "Failed to find target = '%s' for test 7\n", target);
        exit(1);
      }
      break;
  }

  x = F.getSequence();
  fprintf(stderr, "%s\n%s\n", x->header(), x->sequence());
  delete x;
}


void
test3(void) {
  FastASequenceInCore  *x;

  FastAFileWrapper  F(filename);

  int i=0;

  while (!F.eof()) {
    fprintf(stderr, "%d\n", i++);
    x = F.getSequence();
    fprintf(stderr, "%s\n", x->header());
  }
}




int
main(int argc, char **argv) {

  if (argc == 1) {
    fprintf(stderr, "usage: %s <file> <test> <test> <...>\n", argv[0]);
    fprintf(stderr, "tests:\n");
    fprintf(stderr, "       1 == print every sequence in <file>, without using an index\n");
    fprintf(stderr, "       2 == print every sequence in <file>, using an index (basic)\n");
    fprintf(stderr, "       3 == print every sequence in <file>, using an index (id's)\n");
    fprintf(stderr, "       4 == print every sequence in <file>, using an index (names)\n");
    fprintf(stderr, "       5 == print the next arg as iid (basic index)\n");
    fprintf(stderr, "       6 == print the next arg as id (id-only index)\n");
    fprintf(stderr, "       7 == print the next arg as id (defline index)\n");
    fprintf(stderr, "       8 == eof test - read, but don't print, until eof\n");
  }

  filename = argv[1];

  for (int i=2; i<argc; i++) {
    int x = atoi(argv[i]);
    switch (x) {
      case 1:
      case 2:
      case 3:
      case 4:
        test1(x);
        break;
      case 5:
      case 6:
      case 7:
        test2(x, argv[++i]);
        break;
      case 8:
        test3();
        break;
    }
  }


#if 0
  F.find(0);
  x = F.getSequence();
  fprintf(stdout, "h=%s\ns=%s\n", x->header(), x->sequence());
  delete x;

  F.find(1);
  x = F.getSequence();
  fprintf(stdout, "h=%s\ns=%s\n", x->header(), x->sequence());
  delete x;

  F.find(0);
  x = F.getSequence();
  fprintf(stdout, "h=%s\ns=%s\n", x->header(), x->sequence());
  delete x;

  x = F.getSequence();
  fprintf(stdout, "h=%s\ns=%s\n", x->header(), x->sequence());
  delete x;
#endif
}





int crash()
{
  char * p = NULL;
  *p = 0;
  return 0;
}

int foo4()
{
  crash();
  return 0;
}

int foo3()
{
  foo4();
  return 0;
}

int foo2()
{
  foo3();
  return 0;
}

int foo1()
{
  foo2();
  return 0;
}

int main(int argc, char ** argv)
{
  struct sigaction sigact;

  sigact.sa_sigaction = crit_err_hdlr;
  sigact.sa_flags     = SA_RESTART | SA_SIGINFO;

  //  Don't especially care if this fails or not.
  sigaction(SIGSEGV, &sigact, NULL);

  foo1();

  exit(EXIT_SUCCESS);
}

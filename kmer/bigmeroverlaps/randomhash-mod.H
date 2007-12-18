
class randomHash {
public:
  randomHash() {
    mpz_init2(m, 1024);
    mpz_init2(c, 1024);
    mpz_init2(p, 1024);
    mpz_init2(q, 1024);
    gmp_randinit_mt(gmprandomstate);
  };
  ~randomHash() {};

  u32bit   _hashSize;

  mpz_t    m;

  mpz_t    c;
  mpz_t    p;
  mpz_t    q;

  char     str[1024];

  gmp_randstate_t  gmprandomstate;

  //  We want to find a prime p that is as close to 2^width, but no bigger.
  //
  void     setHashSize(u32bit width) {
    _hashSize = width;
    resetRandomness();
  };


  void     resetRandomness(void) {


    ////////////////////////////////////////
    //
    //  p <- a prime bigger than the size of the universe -- 4^100!
    //
    //  We let p = 4^128 + random(2^32), then pick the next prime.
    //
    mpz_ui_pow_ui(c, 2, 32);
    mpz_urandomm(c, gmprandomstate, c);
    mpz_ui_pow_ui(p, 4, 128);
    mpz_add(p, p, c);

    u32bit isprimeloops = 1;
    while (isprimeloops > 0) {
      mpz_add_ui(p, p, 1);
      if (mpz_probab_prime_p(p, 256) == 0)
        isprimeloops--;
    }

    //fprintf(stderr, "Found prime p %s\n", mpz_get_str(str, 10, p));



    ////////////////////////////////////////
    //
    //  1 <= c <= p-1
    //
    mpz_sub_ui(c, p, 2);
    mpz_urandomm(c, gmprandomstate, p);
    mpz_add_ui(c, c, 1);

    //fprintf(stderr, "Found multiplier c %s\n", mpz_get_str(str, 10, c));


    ////////////////////////////////////////
    //
    //  q <- 2^width + 1
    mpz_ui_pow_ui(q, 2, _hashSize);
    mpz_add_ui(q, q, 1);

    //  2 if definitely prime, 1 if possibly, 0 if composite
    //
    //  Proven worst case probability of 0.25 that a strong
    //  pseudoprime will be falsely labeled as prime, so we do a few
    //  iterations to lower the probability.
    //
    isprimeloops = 1 + mtRandom32(mt) % 4096;
    while (isprimeloops > 0) {
      mpz_sub_ui(q, q, 2);
      if (mpz_probab_prime_p(q, 256) == 0)
        isprimeloops--;
    }

    //fprintf(stderr, "Found prime q %s\n", mpz_get_str(str, 10, q));
  };

  u64bit   hash(kMer x) {
    //  ((c * x) mod (p)) mod q)

    //  Copy kMer x into m.
    x.setMPZ(m);

    mpz_mul(m, c, m);
    mpz_mod(m, m, p);
    mpz_mod(m, m, q);

#ifdef TRUE64BIT
    return(mpz_get_ui(m) & u64bitMASK(_hashSize));
#else
    u64bit  r2 = mpz_get_ui(m);
    mpz_tdiv_q_2exp(m, m, 32);
    u64bit  r1 = mpz_get_ui(m);
    return(((r1 << 32) | r2) & u64bitMASK(_hashSize));
#endif
  };
};


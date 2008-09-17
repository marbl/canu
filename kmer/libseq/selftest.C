

  {
    seqFile  *SF = openSeqFile(argv[1]);

    fprintf(stdout, "source '%s' of type '%s' has "u32bitFMT" sequences.\n",
            SF->getSourceName(), SF->getFileTypeName(), SF->getNumberOfSequences());

    fprintf(stdout, "getSequenceLength() vs getSequence(full)\n");
    {
      char  *h = 0L;
      char  *s = 0L;
      u32bit hLen=0, hMax=0;
      u32bit sLen=0, sMax=0;

      for (u32bit sid=0; sid<SF->getNumberOfSequences(); sid++) {
        SF->getSequence(sid, h, hLen, hMax, s, sLen, sMax);

        if ((strlen(s) != SF->getSequenceLength(sid)) ||
            (strlen(s) != sLen) ||
            (SF->getSequenceLength(sid) != sLen)) {
          fprintf(stdout, "length differ for sid="u32bitFMT" h='%s' strlen(s)=%d sLen="u32bitFMT" getSequenceLength()="u32bitFMT"\n",
                  sid, h, strlen(s), sLen, SF->getSequenceLength(sid));
        }
      }

      delete [] h;
      delete [] s;
    }


    fprintf(stdout, "getSequenceLength() vs getSequence(part)\n");
    {
      char  *p = new char [128 * 1024 * 1024];

      for (u32bit sid=0; sid<SF->getNumberOfSequences(); sid++) {
        SF->getSequence(sid, 0, SF->getSequenceLength(sid), p);

        if (strlen(p) != SF->getSequenceLength(sid)) {
          fprintf(stdout, "length differ for sid="u32bitFMT" strlen(s)=%d getSequenceLength()="u32bitFMT"\n",
                  sid, strlen(p), SF->getSequenceLength(sid));
        }
      }

      delete [] p;
    }




  return(0);
}


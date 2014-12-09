

  {
    seqFile  *SF = openSeqFile(argv[1]);

    fprintf(stdout, "source '%s' of type '%s' has "uint32FMT" sequences.\n",
            SF->getSourceName(), SF->getFileTypeName(), SF->getNumberOfSequences());

    fprintf(stdout, "getSequenceLength() vs getSequence(full)\n");
    {
      char  *h = 0L;
      char  *s = 0L;
      uint32 hLen=0, hMax=0;
      uint32 sLen=0, sMax=0;

      for (uint32 sid=0; sid<SF->getNumberOfSequences(); sid++) {
        SF->getSequence(sid, h, hLen, hMax, s, sLen, sMax);

        if ((strlen(s) != SF->getSequenceLength(sid)) ||
            (strlen(s) != sLen) ||
            (SF->getSequenceLength(sid) != sLen)) {
          fprintf(stdout, "length differ for sid="uint32FMT" h='%s' strlen(s)=%d sLen="uint32FMT" getSequenceLength()="uint32FMT"\n",
                  sid, h, strlen(s), sLen, SF->getSequenceLength(sid));
        }
      }

      delete [] h;
      delete [] s;
    }


    fprintf(stdout, "getSequenceLength() vs getSequence(part)\n");
    {
      char  *p = new char [128 * 1024 * 1024];

      for (uint32 sid=0; sid<SF->getNumberOfSequences(); sid++) {
        SF->getSequence(sid, 0, SF->getSequenceLength(sid), p);

        if (strlen(p) != SF->getSequenceLength(sid)) {
          fprintf(stdout, "length differ for sid="uint32FMT" strlen(s)=%d getSequenceLength()="uint32FMT"\n",
                  sid, strlen(p), SF->getSequenceLength(sid));
        }
      }

      delete [] p;
    }




  return(0);
}


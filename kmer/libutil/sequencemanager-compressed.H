
CompressedSequenceManager

Stores loaded sequences in a compressed format (bzip2, 2-bit or 3-bin
encoding, depending on sequence length, and presence of N's).  All
sequences in the file are loaded at the start (option for lazy?), and
compressed then.

When a sequence is requested, a FastABuffer is created and filled
with the sequence.

We could do this offline, too -- store the file compressed (same way),
on access, it's decompressed.

I think the sole-user of sequencemanager is sim4db.  SequenceManager
uses FastABuffers, but doesn't return them.  The sequences returned
are owned by the SeqMgr.

SeqMgrCmp stores it's sequences (headers) compressed.  When requested,
it will uncompress them and return a FastABuffer.  The buffer is then
owned by the caller.  It can be recycled, of course.

Prototypes:

bool     getSequence(FastABuffer &B, u32bit i);
u32bit   numberOfSequences(void);

Constructors:

CompressedSequenceManager(char *filename);

Do we want to be lazy or not?  Load everything at construct time,
or wait to load until we need it?  If we load at init time, we
should keep a pre-compressed image on disk.  If we're lazy, we
should compress on the fly (or maybe still load from compressed?)

For ESTs, this is not a very painful decision, especially if we need
them over and over and over again.  Genomic is harder, as we typically
do not need it over and over again, and there is no need to keep it in
core after we're done with it.

So, we have to ask why we need this?  Why?  search doesn't need it;
sim4 might make use of it, but we should check how many ESTs are used
more than once per file.

Bigger fruit to pick:

1) swap the table building sequences in searchGENOME.

2) sim4faster -- table build, merge (greedy and abort/msp-link),
return of result, driver.

3) play with mRNA in ESTmapper; pick nice abort threshol

4) doc ESTmapper

5) think about alternate splice -- database load, reading from xml


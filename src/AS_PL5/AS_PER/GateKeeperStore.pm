package AS::PER::GateKeeperStore;

=head1 Name

  AS::PER - Perl binding to Celera gate keeper store.

=head1 SYNOPSIS

  use AS::PER;
  my $store = AS::PER::GateKeeperStore->open("./mystore.gkpStore");
  my $iid = $store->UIDtoIID(26);
  my $frag = $store->getFrag($iid, $AS::PER::FRAG_SEQ);
  my $mate_iid = $frag->{"MateIID"};
  my $seq = $frag->{"Sequence"};


=head1 DESCRIPTION

This module is a binding to the Celera gate keeper store. It supports
only reading, not writing, a store.

=head1 METHODS

=over

=item C<< AS::PER::GateKeeperStore->open >>

Open a store. It takes one argument, the path to the store directory.

=item Methods for a store

  NumLibraries()
  NumFragments()
  getFrag(iid, flags=0)
  getLibrary(iid)
  UIDtoIID(uid)
  UIDtoIID_type(uid)

The C<flags> argument to C<getFrag> should be an OR combinaison of the
following flags:

  AS::PER::FRAG_SEQ
  AS::PER::FRAG_QLT
  AS::PER::FRAG_HPS
  AS::PER::FRAG_SRC
  AS::PER::FRAG_ALL

C<UIDtoIID> returns only the value of the IID. C<UIDtoIID_type>
returns an array containing the IID and the type. The type is one of
the following constant:

  AS::PER::IID_UNK
  AS::PER::IID_BAT
  AS::PER::IID_FRG
  AS::PER::IID_LIB

A text description can be obtain from the hash C<%AS::PER::IIDType>
which maps the previous constant to strings ("UNK", "BAT", etc.).

=item Fragment

A fragment behaves like a hash with the following keys available:

  UID
  IID
  MateIID
  LibraryIID
  IsDeleted
  IsNonRandom
  SequenceLength
  QualityLength
  HPSLength
  SourceLength
  Sequence
  Quality
  HPS
  Source

Note that the keys C<Sequence>, C<Quality>, C<HPS> and C<Source> are
valid only if the corresponding flag was passed to the C<getFrag>
method used to retreive the fragment.

=item Library

A library behaves like a hash with the following keys available:

  comment
  mean
  stddev

=back

=head1 LICENSE

This is released under the GPL License.

=head1 AUTHOR

Guillaume Marçais <gus@math.umd.edu>

=head1 SEE ALSO

L<http://wgs-assembler.sourceforge.net/>

=cut


# Extra code for Parser
package AS::MSG::Parser;
use Scalar::Util qw(blessed);

sub open {
    my ($pkg, $h) = @_;
    my $self = AS::MSG::Parser->new();
    if(defined($self)) {
        $self->initialize();
        $self->conf($h);
    }
    return $self;
}

sub conf {
    my ($self, $h) = @_;

    $self->swig_in_set($h->{"IN"}) if(defined($h->{"IN"}));
    $self->swig_out_set($h->{"OUT"}) if(defined($h->{"OUT"}));
    $self->swig_include_set($h->{"INCLUDE"}) if(defined($h->{"INCLUDE"}));
    $self->swig_exclude_set($h->{"EXCLUDE"}) if(defined($h->{"EXCLUDE"}));
    $self->swig_filter_mode_set(1) if($h->{"FILTER_MODE"});
}

sub next_message2 {
    my $self = shift;
    my $m = $self->next_message();
    my %h;

    tie(%h, AS::MSG::MesgWrapper, $m, []);
    return \%h;
}

# Wrapper class around a message to ease memory management issu.  It
# keeps a reference to every object stored in the message, to ensure
# they won't disappear before the message is written or destroyed.
package AS::MSG::MesgWrapper;
sub TIEHASH {
    my ($class, $m, $refs) = @_;
    my $self = { MESG => $m,
                 REFS => $refs,
             };
    return bless $self, $class;
}

sub CLEAR {}

sub EXISTS {
    my ($self, $key) = @_;
    return exists $self->{MESG}->{$key};
}

sub FIRSTKEY {
    my ($self) = @_;
    my $a = keys %{$self->{MESG}};

    each %{$self->{MESG}};
}

sub NEXTKEY {
    my ($self, $key) = @_;

    each %{$self->{MESG}};
}

sub FETCH {
    my ($self, $field) = @_;
    my $res = $self->{MESG}{$field};

    if(ref($res) eq 'ARRAY' && $#$res >= 0 && blessed($res->[0]) &&
       $res->[0]->isa(AS::MSG)) {
        for(my $i = 0; $i <= $#$res; $i++) {
            my %h;
            tie(%h, AS::MSG::MesgWrapper, $res->[$i], $self->{REFS});
            $res->[$i] = \%h;
        }
    }

    return $res;
}

sub STORE {
    my ($self, $field, $value) = @_;
    my $cmd = "swig_" . $field . "_set";
    my (@ary, $res, $tmp);

    if(ref($value) eq 'ARRAY') {
        foreach my $obj (@{$value}) {
            $tmp = tied(%$obj);
            if(ref($obj) eq 'HASH' && ($tmp = tied(%$obj)) &&
               ::blessed($tmp) && $tmp->isa(AS::MSG::MesgWrapper)) {
                push(@ary, $tmp->message());
            } else {
                push(@ary, $obj);
            }
        }
        $res = $self->{MESG}->$cmd(\@ary);
    } else {
        $res = $self->{MESG}->$cmd($value);
    }

    if(defined($res)) {
        push(@{$self->{REFS}}, $res);
    }
}

sub clear_refs {
    my $self = shift;

    $self->{REFS} = ();
}

sub message {
    my $self = shift;

    return $self->{MESG};
}

1;

__END__

=head1 NAME

AS::MSG - Perl binding to Celera message file parser.

=head1 SYNOPSIS

    use AS::MSG;
    my $parser = AS::MSG::Parser->open();
    my $mesg = $parser->next_message();
    print($parser->message_type(), "\n");

=head1 DESCRIPTION

This module is a binding to the Celera message file parser. I.e., it
can parse the files used by Celera (like the .frg, .cgi, .asm, etc.),
using the parser written in C that Celera itself uses.

A message file contains a list of messages, each message containing
some fields and eventually other embedded messages.

Every time the method C<next_message()> is called, a complete message,
including all its embedded messages, is parsed and stored in the
cache. An object pointing to this cache is returned. Then, the fields
of the message can be queried similarly to a hash. For example:

  my $mesg = $parser->next_message();
  $mesg->{"acc"}      # Get the accession number of the message

=head1 Methods

The following methods are defined:

=over

=item C<< AS::MSG::Parser->open >>

Returns a new parser. It takes no arguments or a hash of options. See
L<Options> for the list of options. By default, the parser will read
from standard input and write to standard output.

=item C<< $parser->next_message() >>

Parse the next message from the stream and store it in the cache. It
returns an object or C<undef> if the end of file was reached. To read
all the messages in a file, one can do:

  while(my $mesg = $parser->next_message()) {
    # Do something with $mesg
  }

=item C<< $parser->message_type() >>

Returns a 3 letter string giving the type of the message in the cache
(E.g. AFG, IUM, etc.)

=item C<< $parser->write_message() >>

Write the current cached message to the out stream.

=back

=head1 Options

=over

=item IN

Specify the input stream to read from. It defaults to stdin. E.g.

  open(my $io, "<", "file.asm");
  my $parser = AS::MSG::Parser->open({IN => $io});

=item OUT

Specify the output stream to write to. It defaults to stdout. E.g.

  open(my $io, ">", "mod_file.asm");
  my $parser = AS::MSG::Parser->open({OUT => $io});

=item INCLUDE

List of message types to return. If the list is empty, then
C<next_message()> will return after each message. If the list is not
empty, then C<next_message()> will return only after a message whose
type is in the list (or end of file).

The size of list of types is limited to 16. E.g.:

  my $parser = AS::MSG::Parser->open({INCLUDE => ["MDI", "AFG"]});

=item EXCLUDE

List of message types to not return. C<next_message()> will return
after a message if its type is not in this list, or this list is
empty. The size of the list is limited to 16.

  my $parser = AS::MSG::Parser->open({EXCLUDE => ["UPS", "SCF"]});

=item FILTER_MODE

Act as a filter (think UNIX filter). I.e., every messages which does
not match the INCLUDE list or the EXCLUDE list is written to the
output stream automatically.

E.g.: modify the accession number of all AFG messages, and exclude all MDI messages, copying from stdin to stdout.

  my $parser = AS::MSG::Parser->open({INCLUDE => ["AFG"],
                                      EXCLUDE => ["MDI"],
                                      FILTER_MODE => 1});
  while(my $m = $parser->next_message()) {
    $m->{"acc1"} = 0;
    $parser->write_message();
  }

=back

=head1 Sub-messages

The sub-messages are returned as an array of messages. For example, if
the message is of type IUM, it has IMP sub-messages. Let print the
number of them and the "mid" field of each:

  my $mesg = $parser->next_message();
  my $imps = $mesg->{"IMPs"};
  print("Nb imps: ", scalar(@{$imps}), "\n");
  foreach my $imp (@{$imps}) {
    print("mid: ", $imp->{"mid"});
  }

=head1 Field names


All fields in a message can be retrieved from the cache like a hash
member. For example, if "act" is a field, one can retrieve its value
with

  my $m = $parser->next_message();
  my $act = $m->{"act"};

=head2 General rule

The fields are available via their name in the C structure or the 3
letter acronym as used in the message file. There is an almost
one-to-one correspondance between the 3 letter acronym and the C
structure field name (see Exceptions below). The C header file with
all those structure is src/AS_MSG/AS_MSG_pmesg.h.

For example, a FRG message has a field "version" wich is also
accessible as "ver", the acronym used in the message file. I.e.

  $parser->{"version"} == $parser->{"ver"}

=head2 Sub-messages

The sub-messages of a message are not preceded by a field name. They
can be accessed via the C structure field name (which does not follow
any regular naming scheme) or, by convention, the type of the
sub-message (capitalized) followed by an "s". E.g., "IMPs" can be used
to access the field "f_list", which is the list of sub-message of type
IMP in a IUM message.

=head2 Exceptions

Some field in a message file correspond to 2 fields in the C
structure. This is the case for the "acc:" field of some message which
contains 2 values, written in "eaccession" and "iaccession" in the C
structure. In this case, they are also aliased as "acc1" and
"acc2". The messages for which it happens is: IFG, AFG, UTG, CCO, SCF.

Similarly, in an MDI messages, the "ref:" field is split into
"erefines" and "irefines", and they are aliases as "ref1" and "ref2".

=head2 Example

A complete example, the SCF message, which C structure is given by:

  typedef struct {
    CDS_UID_t           eaccession;
    IntScaffold_ID      iaccession;
    int32               num_contig_pairs;
    SnapContigPairs     *contig_pairs;
  } SnapScaffoldMesg;

And would appear in a file as:

  {SCF
  acc:(1,2)
  noc:2
  {SCP
  ...
  }
  {SCP
  ...
  }
  }

Here, both the "eaccession" and "iaccession" fields are written under
one field in the message file (C<acc:>). By convention, "acc1" and
"acc2" will be aliases for "eaccession" and "iaccession"
respectively.

Then, "num_contig_pairs" field can be accessed as "noc" (general rule).

And the sub-messages of type SCP can be accessed with "contig_pairs" or "SCPs"

=head1 Introspection

The method C<members> of a message returns a hash reference with the
association between the 3/4 letter field name to the C struct field
name. For example, if the message C<$m> is an IUM (IntUnitigMesg
class), the following code:

  my $h = $m->members;
  foreach my $k (keys %{$h}) {
    printf("$k\t%s\n", $h->{$k});
  }

Returns:
  len     length
  cns     consensus
  for     forced
  src     source
  nfr     num_frags
  sta     status
  acc     iaccession
  cov     coverage_stat
  IMPs    f_list
  qlt     quality

The messsage object also responds to C<keys> or C<each>, like a
hash. So:

  print(join(" ", keys %{$m}), "\n");

or

  while(my ($key, $val) = each %{$m}) {
    print($key, ' -> ', $val, "\n");
  }

will work as expected.

=head1 Modifying messages

=head2 The easy way

The easy way to modify messages (but with a performance penality) is
to read them with C<next_message2()>. Then, the standard hash syntax
applies. For example, the following code for each IUM message changes
the coverage to 0, changes the source to "Bogus source", shortens the
list of IMP submessages to 1 element, with a delta list of C<[1, 2, 3]>.

   while(my $m = $parser->next_message2()) {
     $m->{"cov"} = 0;
     $m->{"src"} = "Bogus source";
     my $imps = $m->{"IMPs"};
     splice(@{$imps}, 1);
     $imps->[0]{"delta"} = [1, 2, 3];
     $m->{"IMPs"} = $imps;

     $parser->write_message();
   }

If only reading is performed, the normal C<next_message()> is
preferred to avoid the performance hit (which occurs even for reading).

=head1 The hard way

For some odd memory management issues, modifying messages is a little
more tricky. The C<next_message2()> method wraps the messages in a
class that hides the memory management issues discussed
henceforth. From this point, we assume that the message C<$m> was
obtained with C<next_message()>.

For numerical fields (integer or float), there is no problem, it
suffices to write to the message with the familiar hash syntax. E.g.

  $m->{"acc"} = 0;

For strings or arrays, the hash syntax must NOT be used, but rather
the explicit setter method. The setter for a field called "acc" will
be "swig_acc_set".

Moreover, the return value of the setter must be stored until the
message is written (for example, store it in an array).

For example, the following code is equivalent to the one in the
previous section, but not nearly as pretty.

  my @refs;
  while(my $m = $parser->next_message()) {
    @refs = ();
    $m->{"cov"} = 0;
    push(@refs, $m->swig_src_set("Bogus source"));
    my $imps = $m->{"IMPs"};
    splice(@{$imps}, 1);
    push(@refs, $imps->[0]->swig_delta_set([1, 2, 3]));
    push(@refs, $m->swig_IMPs_set($imps));

    $parser->write_message();
  }

=head2 Caveat

For arrays, the setter method updates automatically the length field
associated with it (for example, in the previous code snippet,
C<swig_delta_set> modified the C<del> field as well as
C<dln>). Sometimes, many fields containing arrays share the same
length fields. The setter does not check that all the lengths
agree. The programmer is responsible to enter coherent data.

In some instances, a string also has a length associated with it. For
example, an IUM message has a C<cns> and C<qlt> fields (for consensus
sequence and quality scores), and a C<len> field giving the length of
the sequence. The setter function for C<cns> and C<qlt> will not
modify the C<len> field. This again is the responsability of the
programmer.

=head1 General caveat

The parser does not play well with Perl IO. I.e., if the parser write
messages to stdout, writting to stdout with the Perl C<print> is not
advised. The results are undetermined, and most likely some
information will be lost.

Keep in mind that the message returned by C<next_message()> referes to
the cached message. It becomes invalid after another call to
C<next_message()>. For example, the following is likely to give strange
results (and may crash).

  my $m1 = $parser->next_message();
  my $m2 = $parser->next_message();
  print($m1->{"cns"}, "\n");

When C<$m1> is used, it may not refer to a valid zone in memory and
hell could ensue.

=head1 LICENSE

This is released under the GPL License.

=head1 AUTHOR

Guillaume Marçais <gus@math.umd.edu>

=head1 SEE ALSO

L<http://wgs-assembler.sourceforge.net/>

=cut

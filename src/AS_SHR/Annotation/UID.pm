# $Id: UID.pm,v 1.3 2008-06-27 06:29:21 brianwalenz Exp $

# Copyright (c) 2003, The Center for the Advancement of Genomics. All rights reserved.

=head1 NAME

UniqueId.pm - A utility for managing unique id's

=head1 VERSION

This document refers to version 1.01 of UniqueId.pm, released 04, 29, 2003.

=head1 SYNOPSIS

=head1 DESCRIPTION

=head2 Overview

Wraps the interface to a uid server.

=cut

package Annotation::UID;

use strict;
use Log::Log4perl;
use fields qw(
  batchSize
  current
  top
  namespace
);

my $log;

=head1 METHODS

=cut

=head2 <Annotation::Sequence> new()

Constructor; Create a new instance of Sequence.

=cut

sub new($$;$) {
  my ($that,$size,$namespace) = @_;

  $log = Log::Log4perl->get_logger($that) unless defined($log);
  my Annotation::UID $this = fields::new( $that );
  # initialize default values
  $size = 100 unless defined($size);
  $size =~ /(\d+)/;
  my $num = $1;
  $log->logconfess("UID batch size must be a positive whole number ($size is invalid)") if $num ne $size;
  $this->{batchSize} = $size;
  $this->{current} = undef;
  $this->{top} = undef;
  $this->{namespace} = $namespace if defined $namespace;

  return $this;
}

=head1 METHODS

=cut

=head2 <long> incrUID()

Return a new uid.

=cut

sub incrUID($) {
  (my Annotation::UID $this) = @_;

  $log = Log::Log4perl->get_logger($this) unless defined($log);

  if (defined($this->{current}) && $this->{current} != $this->{top}) {
    $this->{current}++;
  } else {
    # need a new batch of UIDs
    my $url = "http://guid.jtc.jcvsf.org:8080/guid/GuidClientServer?Request=GET&Size=$this->{batchSize}";
    $url .= "&Namespace=$this->{namespace}" if defined $this->{namespace};
    my $result = `wget -q -O - \"$url\"`;
    $result =~ s/\s+//gs;
    $log->debug("Result from UID server:\n$result");
    if ($result =~ /SUCCESS/ && $result =~ />(\d+)</s) {
      $this->{current} = $1;
      $this->{top} = $1 + $this->{batchSize} - 1;
    } else {
      $log->logconfess("Error retrieving UIDs (\"$result\")");
    }
  }

  return $this->{current};
}

=head2 <long> getUID()

Gets the current (most recent) UID.

=cut

sub getUID($) {
  my ($this) = @_;

  return $this->{current};
}


=head2 <long> getBatchSize()

Gets the batch size for this UID server object.

=cut

sub getBatchSize($) {
  my ($this) = @_;

  return $this->{batchSize};
}

1;

__END__

=head1 ENVIRONMENT

Not applicable.

=head1 DIAGNOSTICS

Uses standard debugging and logging tools.

=head1 BUGS

There is no code here to access an external UID server.

=head1 SEE ALSO

Log::Agent

=head1 AUTHOR(S)

  Doug Rusch
  The Center for the Advancement of Genomics
  1901 Research Blvd, 6th Floor
  Rockville, MD 20850

=head1 COPYRIGHT

Copyright (c) 2003, The Center for the Advancement of Genomics. All Rights Reserved.


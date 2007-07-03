#! /usr/local/bin/perl -w

# Copyright (c) 2004  The Institute for Genomic Research.  All rights reserved.

use strict;
use POSIX;
use TIGR::Foundation;
use File::Basename;
use Config::IniFiles;

# This program creates an instance configuration as the union / intersection
# of three input files: the global configuration, the pipeline configuration,
# and the user supplied configuration.  This configuration file is written
# to a temporary location.

# Define global variables.
my $PROGNAME = "/usr/local/common/crs/common/execPipeline";
my $GLOBAL_CONFIG = "/usr/local/common/crs/etc/global.conf";
my $PIPELINE_CONFIG = "/usr/local/common/crs/etc/buildTiling.conf";

# Set a PERLLIB variable as necessary.
if ( ! defined ( $ENV{'PERLLIB'} ) ) {
   $ENV{'PERLLIB'} = "/usr/local/common/crs/lib";
}
else {
   $ENV{'PERLLIB'} = $ENV{'PERLLIB'} . ":" . "/usr/local/common/crs/lib";
}

my $tf_obj = new TIGR::Foundation;              # init TF
# Load the TIGR Foundation information from the programs themselves.
my $helptext = `$PROGNAME -help 2>&1`;
my $versionstring = `$PROGNAME -version 2>&1`;
my @depends = split /\n/, `$PROGNAME -depend 2>&1`;
# Doctor the help string.
$helptext =~ s/execPipeline/buildTiling/g;
# Doctor the version string.
$versionstring =~ s/^\w+ //;
chomp $versionstring;
$tf_obj->setVersionInfo($versionstring);
$tf_obj->setHelpInfo($helptext);
$tf_obj->addDependInfo(@depends);

# Create a temporary file.
my $tmpdir = $ENV{'ACTMPDIR'};
if ( ! defined ( $tmpdir ) ) {
   $tmpdir = "/local/aserver/AC_WORK";
}
my @new_argv = @ARGV;                           # receive ARGV
my $tmpfile = $tmpdir . "/" . basename(POSIX::tmpnam());
my $prefix = undef;
$tf_obj->TIGR_GetOptions( "l=s" => \$prefix );
# Check the size of ARGV after TF processing.  If it's got arguments,
# pop off the last argument for use as a configuration file.
my $user_config;
if ( scalar(@ARGV) > 0 ) {
   $user_config = pop @ARGV;
   # Pop off the ARGV from the stored version.  It will be at the end.
   pop @new_argv;
}

# Start with the global configuration, then the pipeline config,
# then the user supplied config as applicable. Unify the configurations.
my %Pipeline = ();
my @cf_files = ( $GLOBAL_CONFIG, $PIPELINE_CONFIG );
if ( $user_config ) {
   push @cf_files, $user_config;
}
foreach my $cf ( @cf_files ) {
   my $conf_obj = new Config::IniFiles ( -file => $cf );
   if ( ! $conf_obj ) { 
      die "Failed to load configuration file \'$cf\'.\n";
   }
   my @sections = $conf_obj->Sections();
   foreach my $sect_name ( @sections ) {
      if ( ! defined ( $Pipeline{$sect_name} ) ) {
         $Pipeline{$sect_name} = {};         # Create a new section.
      }
      foreach my $parm ( $conf_obj->Parameters($sect_name) ) {
         $Pipeline{$sect_name}->{$parm} = $conf_obj->val($sect_name, $parm);
      }
   }
}
# Write the new configuration to the temporary file.
if ( ! open ( TMPFILE, ">" . $tmpfile ) ) {
   die "Failed to open temporary configuration file \'$tmpfile\' for " .
      "writing.\n";
}
foreach my $sect_name ( keys %Pipeline ) {
   print TMPFILE "[", $sect_name, "]", "\n";
   foreach my $parm ( keys %{$Pipeline{$sect_name}} ) {
      print TMPFILE $parm, "=", $Pipeline{$sect_name}->{$parm}, "\n";
   }
}
if ( ! close (TMPFILE) ) {
   die "Failed to close temporary configuration file \'$tmpfile\'.\n";
}
# Add new configuration file to invocation.
push @new_argv, $tmpfile;

exit exec "$PROGNAME @new_argv";

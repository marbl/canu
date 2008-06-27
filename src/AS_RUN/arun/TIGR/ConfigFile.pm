# $Id: ConfigFile.pm,v 1.2 2008-06-27 06:29:19 brianwalenz Exp $

# Copyright @ 2002 - 2010 The Institute for Genomic Research (TIGR).
# All rights reserved.
#
# This software is provided "AS IS".  TIGR makes no warranties, express or
# implied, including no representation or warranty with respect to the
# performance of the software and derivatives or their safety,
# effectiveness, or commercial viability.  TIGR does not warrant the
# merchantability or fitness of the software and derivatives for any
# particular purpose, or that they may be exploited without infringing the
# copyrights, patent rights or property rights of others.
#
# This software program may not be sold, leased, transferred, exported or
# otherwise disclaimed to anyone, in whole or in part, without the prior
# written consent of TIGR.

package TIGR::ConfigFile;
{

=head1 NAME

TIGR::ConfigFile - The TIGR::ConfigFile module provides access to variables
stored in a file that has INI file format

=head1 SYNOPSIS

   use TIGR::ConfigFile;
   my $config_obj = new TIGR::ConfigFile($config_file, $foundation_obj);
   my $value = $config_obj->getOption($option, $section);

=head1 DESCRIPTION

This module provides an object oriented interface for retrieving
configuration information.

=cut

   BEGIN {
      require 5.006_00;
   }

   use strict;
   use TIGR::Foundation;

   our $REVISION = (qw$Revision: 1.2 $)[-1];
   our $VERSION = '1.01';
   our $VERSION_STRING = "$VERSION (Build $REVISION)";
   our @DEPEND = ('TIGR::Foundation');

   # debugging scheme
   #
   #   Debugging via the TIGR Foundation uses increasing log levels based on
   #   nesting.  'MAIN' starts at level 1.  Every nest increments the level by
   #   1. Subroutines always start nesting at level 2.  As debugging levels
   #   increase, logging is more verbose.  This makes sense as you log at
   #   greater depth (ie. deeper branching).
   #
   #   The following definitions help emphasize the debugging in the program.
   #

   my $DEBUG_LEVEL_1 = 1;
   my $DEBUG_LEVEL_2 = 2;
   my $DEBUG_LEVEL_3 = 3;
   my $DEBUG_LEVEL_4 = 4;
   my $DEBUG_LEVEL_5 = 5;
   my $DEBUG_LEVEL_6 = 6;
   my $DEBUG_LEVEL_7 = 7;
   my $DEBUG_LEVEL_8 = 8;
   my $DEBUG_LEVEL_9 = 9;

   # default common section
   my $DEFAULT_SECTION_NAME = ".main";

   # the list of regular expressions used to parse the config file using INI
   # file format.
   my $REGEX_EMPTY_STRING_WITH_WHITESPACE = '^\s*$';
   my $REGEX_EQUALS_SIGN = '\=';
   my $REGEX_COMMENT = '^[;#]';
   my $REGEX_SECTION = '\[([[:print:]]+)\]';
   my $REGEX_SUBSTITUTE_PATTERN = '{{([[:print:]]+?)}}';

   # The name for the default common section
   my $MAIN_SECTION = ".main";

   ## prototypes

   sub new($;$$);
   sub rereadConfig();
   sub getOptions(;$);
   sub getSections();
   sub hasSection($);
   sub getOptionNames(;$);
   sub hasOption($;$);
   sub getOption($;$$);
   sub hasConfigFile();
   sub __readOptions();
   sub __substitute($$$);
   sub __substitute_option($$);
   sub __errorHandler($);
   sub __logHandler($$);
   ## implementation documentation

=over

=item $obj_instance = new TIGR::ConfigFile(config_file, foundation_obj,
$error_array_ref)

This method returns a new instance of the TIGR::ConfigFile object. It takes in
a config file name and two optional parameters - a TIGR::Foundation object and
a reference to an array for storing error messages. The config file is
parsed for correctness in this method. If the config file does not parse or
if errors occur during object creation, this method returns undef. Errors in
parsing are written to the error array passed in and the .error file.

=cut

   sub new($;$$) {
      # The class name
      my $pkg = shift;
      # The array of arguments passed to new
      my @method_args = @_;
      # A variable that tracks the errors in this method
      my $error_condition = 1;
      # A reference to a TIGR::ConfigFile object
      my $self = {};
      # telling $self that it is now a TIGR::ConfigFile object
      bless $self, $pkg;
      # a hash reference that stores all the configuration information
      $self->{__all_options} = {};

      # storing the config file name in the object
      if ( scalar (@method_args) > 0 ) {
         my $filename =  shift @method_args;

         if(defined $filename) {
	    $self->{config_file} = $filename;
            $self->__logHandler("Got config file in new()", $DEBUG_LEVEL_4);
	 }
         else {
            $self->__errorHandler("Undefined file name passed in new");
	 }
      }
      else { # an error since the config file name is required
         $self->{config_file} = undef;
         $self->__errorHandler("No config file in new");
         $error_condition = undef;
      }

      # storing the foundation object in the ConfigFile object
      if ( ( scalar (@method_args) > 0 ) &&
           ( ( ref ($method_args[0]) ) =~ /foundation/i ) ) {
         $self->{foundation} = shift @method_args;

         $self->__logHandler("Got TIGR::Foundation in new()", $DEBUG_LEVEL_3);
      }
      else { # the foundation object was not passed, but no error since they
             # it is optional
         $self->{foundation} = undef;
         $self->__logHandler("No TIGR::Foundation in new()", $DEBUG_LEVEL_3);
      }

      if ( ( scalar (@method_args) > 0 ) &&
           ( ( ref ($method_args[0]) ) =~ /array/i ) ) {
         $self->{error_ref} = shift @method_args;
         $self->__logHandler("Got Error ARRAY in new()", $DEBUG_LEVEL_3);
      }
      else {
         $self->{error_ref} = undef;
         $self->__logHandler("No Error ARRAY in new()", $DEBUG_LEVEL_3);
      }

      # check for invocation errors
      if ( ( scalar (@method_args) > 0 ) ) { # too many parameters passed
         $self->__errorHandler("Too many parameters passed to new() ".
                                 "method");
         $error_condition = undef;
      }

      # Reading in the information in the Config file
      if((defined $error_condition) && (!defined ($self->__readOptions()))) {
         $self->__errorHandler("The config file could not be read or ".
                   "parsed");
         $error_condition = undef;
      }
      # returning a reference to a TIGR::ConfigFile object if there were
      # no errors
      return ($error_condition) ? $self : undef;
   }

=item $obj_instance->rereadConfig()

This method purges the options currently buffered, and re-reads options from
the configuration file. It returns 1 on success and undef on failure

=cut

   sub rereadConfig() {
      my $self = shift;
      # A variable that tracks the errors in this method
      my $error_condition = 1;
      # re-initializing the hash that stores the configuration information
      $self->{__all_options} = {};

      # rereading the configuration information from the config file
      if(!defined($self->__readOptions())) {
         $self->__errorHandler("The config file could not be ".
                                 "re-read or parsed");
         $error_condition = undef;
      }
      return ($error_condition) ? 1 : undef;
   }

=item $obj_instance->getOptions(section)

This method takes in an optional parameter section that specifies a section
name in the config file. If the section is found in the config file the method
returns a hash reference of the option value pairs in that section. If no
section name is passed, a hash reference of the option value pairs in the
common section is returned.  The reused options names are substituted in the
returned hash to obtain the actual values of the options in any section. If
the section is not found or there are other errors the function returns undef.

=cut

   sub getOptions(;$) {
      # a reference to a TIGR::ConfigFile object
      my $self = shift;
      # the section name
      my $section = shift;
      # A variable that tracks the errors in this method
      my $error_condition = 1;
      # the hash reference to the stored parameters
      my $ret_hash = undef;

      # if the section name is not specified, use the main section
      if(!defined $section) {
         $section = $MAIN_SECTION;
      }

      if (defined $section) { # when the section name is defined
         if($self->hasSection($section)) { # the section is present in the
                                           # config file
            $ret_hash = $self->{__all_options}->{$section};
         }
         else { # section name not present hence error
            $self->__errorHandler("The section name $section is not ".
                      "present in the config file");
            $error_condition = undef;
         }
      }

      return ($error_condition) ? $ret_hash : undef;
   }

=item $obj_instance->getSections()

This method returns a reference to an array of sections in the config file.
The common section name is called ".main". If there are no sections in the
config file or if there are other errors this method returns undef

=cut

   sub getSections() {
      my $self = shift;
      # the array that holds the section names
      my @section_arr = ();
      # A variable that tracks the errors in this method
      my $error_condition = 1;
      # The return value from this function
      my $return_val = undef;
      # a reference to the hash of all the section names
      my $option_ref = {};

      if(defined($option_ref = $self->{__all_options})) { # if the hash of
                                                          # options is defined
         @section_arr = keys(%$option_ref); # get the list of section names
         $return_val = \@section_arr;
      }
      else {
         $self->__errorHandler("The configuration hash is not defined");
         $error_condition = undef;
      }

      return ($error_condition) ? $return_val : undef;
   }

=item $obj_instance->hasSection(section)

This method takes in a section and returns 1 if that section is present in
the config file and undef otherwise. The common section name is called ".main"

=cut

   sub hasSection($) {
      my $self = shift;
      # the section name
      my $section = shift;
      # A variable that tracks the errors in this method
      my $error_condition = 1;

      if(defined $section) { # if the section name was passed
         # check if the section name is present in the config file
         if(exists($self->{__all_options}->{$section})) {
            $self->__logHandler("The section $section is present ".
                      "in the config file", $DEBUG_LEVEL_4);
         }
         else { # if the section is not in the config file
            $self->__errorHandler("The section $section is not ".
                      "present in the config file");
            $error_condition = undef;
         }
      }
      else { # if the section name is not defined
         $self->__errorHandler("No section name passed to hasSection");
         $error_condition = undef;
      }
      return ($error_condition) ? 1 : undef;
   }

=item $obj_instance->getOptionNames(section)

This method take in an optional section and returns a reference to an array of
option names in that section. If section is not specified, a reference to an
array of option names in the common section is returned. The common section
can also be represented by setting section to ".main". If the section is not
defined or there are errors the method returns undef

=cut

   sub getOptionNames(;$) {
      my $self = shift;
      # the section name
      my $section = shift;
      # A variable that tracks the errors in this method
      my $error_condition = 1;
      # the return value from this function
      my $return_val = undef;
      # the array of options in this section
      my @option_arr = ();
      # the hash containing the option value pairs for the section passed in
      my $section_ref = {};
      # if the section name is not specified, hence return the options from
      # the common section.
      if(!defined $section) {
         $section = $MAIN_SECTION;
      }

      if(defined $section) { # if the section name is defined
         if(defined($section_ref = $self->{__all_options}->{$section})) {
            @option_arr = keys(%$section_ref ); # assign the array of options
            $return_val = \@option_arr;
         }
         else {
            $self->__errorHandler("Could not obtain the options for section".
                      "$section");
            $error_condition = undef;
         }
      }

      return ($error_condition) ? $return_val : undef;
   }

=item $obj_instance->hasOption(option, section)

This method takes an option and an optional section and returns 1 if the
specified option exists for the specified section. If the section is not
specified, the method checks if the specified option exists in the common
section. The common section can also be represented by setting section to
".main". If the section does not exist or the option does not exist for the
specified section the method returns undef

=cut

   sub hasOption($;$) {
      my $self = shift;
      # The option name for the section
      my $option = shift;
      # The section name specified
      my $section = shift;
      # A variable that tracks the errors in this method
      my $error_condition = 1;

      # if the section name is not specified, use the main section
      if(!defined $section) {
         $section = $MAIN_SECTION;
      }

      if(!defined $option) { # if the option is defined
         $self->__errorHandler("No option name passed to hasOption");
         $error_condition = undef;
      }
      elsif(defined $section) { # if the section is defined
         # check if the section name is present in the config file
         if(defined($self->{__all_options}->{$section})) {
            $self->__logHandler("The section $section is present ".
                      "in the config file", $DEBUG_LEVEL_4);
         }
         else { # if the section is not in the config file
            $self->__errorHandler("The section $section is not ".
                      "present in the config file");
            $error_condition = undef;
         }
         # check if the option exists in the specified section
         if((defined $error_condition) &&
            (!exists($self->{__all_options}->{$section}->{$option}))) {
            $self->__errorHandler("The option $option does not ".
                      "exist in section $section");
            $error_condition = undef;
         }
      }

      return ($error_condition) ? 1 : undef;
   }

=item $obj_instance->getOption(option, section, default)

This method takes in an option name and two optional parameters - a section
name and a default value for the option. If the section is specified, it
returns the value of the option in that section. If the section is not
specified, the value of the option from the common section is returned. The
common section can also be represented by setting section to ".main". The
reused options names are substituted in the returned value to obtain the actual
values of the option in that section. If the option is not found in a section
and the default value is specified, the default value is returned. If the
section is not found or there are other errors, the method returns undef.

=cut

   sub getOption($;$$) {
      my $self = shift;
      # The option name
      my $option = shift;
      # The section name
      my $section = shift;
      # The default value for the option
      my $default = shift;
      # The value of the option in that section
      my $value = undef;
      # A variable that tracks the errors in this method
      my $error_condition = 1;

      if(!defined $option) { # an error since no option name is passed
         $self->__errorHandler("No option name passed to getOption");
         $error_condition = undef;
      }
      else { # if the option is specified
         # if the section name is not specified, use the main section
         if(!defined $section) {
            $section = $MAIN_SECTION;
         }

         # check if the option is defined in the section
         if(!defined($self->hasOption($option, $section))) {
            $self->__logHandler("The option $option is not defined in ".
                      "section $section", $DEBUG_LEVEL_4);
            if(defined $default) { # set to the default value if the option
                                   # is not defined
               $value = $default;
            }
            else { # no default specified
               $self->__logHandler("No default specified for option ".
                         "$option in section $section", $DEBUG_LEVEL_5);
               $error_condition = undef;
	    }
         }
         else { # the option is defined, so take the value of the option in
                # the section
            $self->__logHandler("The option $option is defined in section ".
                      "$section", $DEBUG_LEVEL_4);
            $value = $self->{__all_options}->{$section}->{$option};
         }
      }
      return ($error_condition) ? $value : undef;
   }

=item $obj_instance->hasConfigFile()

This method returns 1 if a config file exists for the object and 0 otherwise.

=cut

   sub hasConfigFile() {
      my $self = shift;
      return $self->{__has_config_file}; # return the flag that specifies if
                                         # the config_file has been set for
                                         # the object
   }

# $obj_instance->__readOptions()

# This method opens the config file and parses the file according to the
# INI file specification. It also reads stores the configuration information
# from the config file in a multilevel hash. The common section of the file is
# called ".main" when storing the options in the common section. The reused
# options names are substituted and stored in the hash to have ready access to
# the actual values of the options. The method returns 1 on success and undef
# otherwise.

   sub __readOptions() {
      my $self = shift;
      # The line number in the config_file
      my $line_num = 0;
      # A variable that tracks the errors in this method
      my $error_condition = 1;
      # The option name in a section
      my $option = undef;
      # The value of an option in a section
      my $value = undef;
      # The array of options
      my @options = ();
      # The array of sections
      my @sections = ();
      # A specific section
      my $section = undef;
      # The config file
      my $config_file = $self->{config_file};
      my %section_hash = ();
      # opening the config file
      if((defined $config_file) && (!open(CONFIGFILE, $config_file))) {
         $self->__errorHandler("Could not open the config_file");
         # setting the config file flag to zero indicating that the file
         # could not be opened.
         $self->{__has_config_file} = 0;
         $error_condition = undef;
      }
      else {
         # The file was opened so set the config file flag to 1
         $self->{__has_config_file} = 1;
      }
      my $line = undef;

      # if a config file has beeen set
      if($self->{__has_config_file}) { # BEGIN_HAS_CONFIG
         my $section = $MAIN_SECTION;
         my $line = undef;

         # Read the config file
         while($line = <CONFIGFILE>) { # BEGIN_READ_CONFIG
	    $line_num = $line_num + 1;
            chomp($line); # remove the new line char from the end

            # Skip, if a blank line.
            if($line =~ /$REGEX_EMPTY_STRING_WITH_WHITESPACE/) {
               next;
            }
            # Skip, if a "comment".
            elsif($line =~ /$REGEX_COMMENT/) {
               next;
            }
            # if the line is a section name
            elsif($line =~ /$REGEX_SECTION/) {
               $section = $1;
               # initialize the hash reference for the section being read
               # check if the section is defined and does not contain only
               # space characters
               if((!defined $section) || ($section =~ /^\s*$/)) {
                  $self->__errorHandler("The section name is not properly ".
                            "specified on line $line_num");
                  $error_condition = undef;
               }
               else { # the section name is defined
                  my $section_check =  $section_hash{$section};
                  if((defined $section_check) && ($section_check == 1)) {
                     $self->__errorHandler("Duplicate section name $section ".
                               "found");
                     $error_condition = undef;
		  }
                  else {
                     $section_hash{$section} = 1;
		  }
                  # initializing the hash reference for the section being read
                  $self->{__all_options}->{$section} = {};
                  $self->__logHandler("Reading the information in ".
                            "the section named $section", $DEBUG_LEVEL_5);
	       }
            }
            # if the line is an option specification
            elsif($line =~ /$REGEX_EQUALS_SIGN/) {
               # obtain the option value pair
               my @option_arr = split(/=/,$line);
               $option = shift @option_arr;
               if((defined $line) && (defined $option)) {
                  $line =~ s/$option=//;
                  $value = $line;
                  if($value eq "") {
		     $self->__logHandler("the value for option $option is ".
                               "empty", $DEBUG_LEVEL_6);
                  }
               }
               # check if both the option and the value are defined and
               # do not contain only space characters
               if((!defined $option) || ($option eq "") ||
                  ($option =~ /^\s*$/) || (!defined $value)){ # a parse error
                  $self->__errorHandler("The option specification is ".
                           "not properly formatted on line $line_num");
                  $error_condition = undef;
               }
               else { # the option specification is correct so store the
                      # information in the hash
                  $option =~ s/\s*//g; # remove the spaces in the option name
                  $value =~ s/^\s*//; # remove the spaces in the beginning of
                                      # the value
                  #check if section is defined
                  if(defined $section) {
                     $section =~ s/\s*//g; # remove the spaces in the section
                                           # name
                     $self->__logHandler("Reading the option $option in ".
                               "section $section" , $DEBUG_LEVEL_6);
                     # store the value for the option in the hash
                     $self->{__all_options}->{$section}->{$option} = $value;
                  }
                  else {# section is not defined
                     $self->__errorHandler("The section name is not defined ".
                               "for option $option");
                     $error_condition = undef;
		  }
               }
            }
            else { # parse error since the line is not valid
                $self->__errorHandler("The config file is ".
                          "not properly formatted on line $line_num");
                $error_condition = undef;
            }
         } # END_READ_CONFIG
      } # END_HAS_CONFIG

      my $option_ref = $self->{__all_options};
      # The number of sections stored
      my $section_num = scalar(keys(%$option_ref));

      # check if there are sections in the config file
      if((!defined $section_num) || ($section_num == 0)) {
         $self->__errorHandler("No sections found in the config file");
         $error_condition = undef;
      }
      # the hash representing the main section
      my $main_options = $self->{__all_options}->{$MAIN_SECTION};
      my @mopts = keys(%$main_options);

      # substituting the reused option names in the common section
      if(!defined($self->__substitute( $MAIN_SECTION , $main_options ,
           $main_options))) {
         $self->__errorHandler("Could not substitute the main ".
                   "options in the config file");
         $error_condition = undef;
      }
      # the array of sections
      @sections = keys(%$option_ref);
      my $section_ref = {};
      my %substitute_options = ();

      # logging the config file information before substitution
      foreach $section (@sections) {
         $self->__logHandler( " Section $section ", $DEBUG_LEVEL_3);
         $section_ref = $self->{__all_options}->{$section};
         @options = keys(%$section_ref);
         foreach $option (@options) {
            $value = $self->{__all_options}->{$section}->{$option};
            if(defined $value) {
               $self->__logHandler( " $option = $value ", $DEBUG_LEVEL_4);
	    }
         }
      }
      $section_ref = {};
      # Perform substitutions for values with text within "{{}}" (curly braces)
      foreach $section (@sections){
         $section_ref = $self->{__all_options}->{$section};

         if ($section eq $MAIN_SECTION) {
	    # Skip -- already handled above.
            $self->__logHandler("skipping main section...", $DEBUG_LEVEL_4);
            next;
         }

         $self->__logHandler( "section $section ",$DEBUG_LEVEL_3);
         my %main_options_hash = %$main_options;
         my %section_ref_hash = %$section_ref;
         # combine the option value pairs from the common section and the
         # section being processed.
         %substitute_options = (%main_options_hash, %section_ref_hash);

         # substituting the reused option names in the section being processed
         if(!defined($self->__substitute( $section, $section_ref,
                               \%substitute_options))) {
            $self->__logHandler("Could not substitute the $section options ".
                      "in the config file", $DEBUG_LEVEL_3);
            $error_condition = undef;
         }
      }

      # logging the config file information after substitution in the
      # log file
      foreach $section (@sections) {
         $self->__logHandler(" Section $section ", $DEBUG_LEVEL_4);
         @options = keys(%$section_ref);
         foreach $option (@options) {
            $value = $self->{__all_options}->{$section}->{$option};
            if(defined $value) {
               $self->__logHandler( " $option = $value ", $DEBUG_LEVEL_5);
	    }
         }
      }

      return ($error_condition) ? 1 : undef;
   }

# $obj_instance->__substitute(section , options , substitute_options)

# This method takes in three parameters - a section name , a reference to a
# hash of options in the section and a reference to a hash of options that are
# used for substitution. This method calls other methods to perform a
# recursive substitution of the reused options names in that section to obtain
# the actual values of the options. It stores the actual values of the options
# in the object hash so they can easily be accessed by getOption. This method
# returns 1 on success and undef otherwise.

   sub __substitute($$$) {
      my $self = shift;
      # The section name
      my $section = shift;
      # a reference to a hash of options in a section
      my $options_ref = shift;
      # a reference to a hash of options that are used for substitution
      my $substitute_options_ref = shift;
      # an option name
      my $option = undef;
      # A variable that tracks the errors in this method
      my $error_condition = 1;
      # the return value from this function
      my $ret_val = undef;

      # iterate through the options in the option hash
      foreach $option (keys(%$options_ref)) {
         # The value of the option
         my $value = $options_ref->{$option};

         $self->__logHandler( "section [$section]  option ".
                   "[$option] value [$value]",  $DEBUG_LEVEL_3);

         # substitute the text within "{{ }}" in the value of the option
         if(!defined($ret_val= $self->__substitute_option(
                                        $value,
                                        $substitute_options_ref))) {
            # the substitution failed hance an error
            $self->__errorHandler( "The option [$option] in section ".
                      "[$section] could not be substituted");
            $error_condition = undef;
         }
         else { # the substitution was successful
            $self->__logHandler("The value of the substituted option ".
                      "[$option] in section [$section] is [$ret_val]",
                      $DEBUG_LEVEL_4);
            # store the substituted value of the option in the object hash
            $self->{__all_options}->{$section}->{$option} = $ret_val;
         }
      }
      return ($error_condition) ? 1 : undef;
   }

# $obj_instance->__substitute_option(text, substitutes)

# This method takes in two parameters - a string val that is to be substituted
# and a reference to a hash of options that are used for substitution. This
# method recursively substitutes parts of the text that lie between "{{ }}".
# This method returns the substituted text on success and undef otherwise.

   sub __substitute_option($$) {
      my $self = shift;
      # The string that is to be sustituted
      my $text = shift;

      # Reference to a hash of options that are used for substitution
      my $substitutes_ref = shift;
      # The part of the text that lies within double curly braces.
      my $text_sub = undef;
      # The value of the part of the text that lies within double curly braces
      my $text_sub_val = undef;
      # A variable that tracks the errors in this method
      my $error_condition = 1;

      if((defined $text) && ($text =~ /$REGEX_SUBSTITUTE_PATTERN/)) {
         # The text contains double curly braces. The part of the text within
         # double curly braces is called replacement text
         $text_sub = $1;

         if(defined $text_sub) {
            $text_sub =~ s/\s*//g;
            $text_sub =~ s/}}$//g;
	    if($text_sub eq "") {
               $self->__errorHandler("The replace text is whitespace");
               $error_condition = undef;
            }

            if(defined ($text_sub_val = $substitutes_ref->{$text_sub})) {
               # the value of the replacement text has been found in the hash
               # of options
               $text =~ s/{{[A-z0-9_]+}}/$text_sub_val/;
            }
            else { # the value of the replacement text was not found in the
                   # hash of options
               $self->__errorHandler("The value of option $text_sub ".
                         "could not be found");
               $error_condition = undef;
            }
         }
      }

      # check if the text has still more double curly braces
      if((defined $text) && ($text =~ /$REGEX_SUBSTITUTE_PATTERN/)) {
         # recursively substitute the partially substituted text
         if(!defined ($text = $self->__substitute_option($text,
                                 $substitutes_ref))) {
            # could not substitute the partially substituted text
            $self->__errorHandler("The call to __substitute_option failed");
            $error_condition = undef;
         }
      }
      return ($error_condition) ? $text : undef;
   }

#$obj_instance->__errorHandler($message);

# This method takes in one parameter - a message to be logged via logError()
# using TIGR::Foundation. The method logs the message via logError. The method
# returns 1 on success and undef if the foundation object is not defined or
# in case of other errors.

   sub __errorHandler($) {
      my $self = shift;
      # getting the message
      my $message = shift;;
      # A variable that tracks the errors in this method
      my $error_condition = 1;

      if ( defined ($message) ) { # if the message was passed

         if ( defined ($self->{foundation}) ) { # if a foundation object is
                                                # defined
            $self->{foundation}->logError($message); # write to the .error
                                                     # file
         }
         else { # foundation object is not defined
            $error_condition = undef;
         }
      }
      elsif(defined ($self->{foundation})) { # incorrect function call
         $self->{foundation}->logError("No message passed to the ".
                                 "errorHandler");
         $error_condition = undef;
      }

      if ( defined ($self->{error_ref}) )  {
         push @{$self->{error_ref}}, $message;
      }

      return ($error_condition) ? 1 : undef;
   }

#$obj_instance->__logHandler($message, $debuglevel);

# This method takes in two parameters - a message to be logged via logLocal()
# using TIGR::Foundation and the debug level for the log message. The method
# logs the message via logLocal. The method returns 1 on success and undef if
# the foundation object is not defined or in case of other errors.

   sub __logHandler($$) {
      my $self = shift;
      # getting the message
      my $message = shift;
      # getting the debug level
      my $level = shift;
      # A variable that tracks the errors in this method
      my $error_condition = 1;

       # if both the message and the debug level were passed
      if ( (defined $message) && (defined $level) ) {
         if ( defined ($self->{foundation}) ) { # if a foundation object is
                                                # defined
            $self->{foundation}->logLocal($message, $level); # write to the log
                                                             # file
         }
         else { # foundation object is not defined
            $error_condition = undef;
         }
      }
      elsif(defined ($self->{foundation})) { # incorrect function call
         $self->{foundation}->logError("Insufficient parameters passed to".
                                 "errorHandler");
         $error_condition = undef;
      }

      return ($error_condition) ? 1 : undef;
   }

=head1 USAGE

To use this module, load the C<TIGR::ConfigFile> package via the
C<use> function. Then, create a new instance of the object via the
C<new()> method, as shown below.

An example script("Test") using this module follows. The C<TIGR::Foundation>
module is included for completeness but does not have to be used.
The configuration file "test_conf" used in this example has an INI format and
contains the following information.

# "test_conf"

    root = /usr/local
    [TEST]
    APPLICATION_ID = TEST
    APPLICATION_NAME = TEST
    APPLICATION_ROOT_DIR = root/pfgrc
    APPLICATION_LOG_SUBDIR = {{APPLICATION_ROOT_DIR}}/log
    APPLICATION_DATA_SUBDIR = {{APPLICATION_ROOT_DIR}}/data
    APPLICATION_LOG_FILESPEC = {{APPLICATION_LOG_SUBDIR}}/test.log

   #!/usr/local/bin/perl -w
   use strict;
   use TIGR::ConfigFile;
   use TIGR::Foundation;

   MAIN:
   {
      # create new foundation object
      my $tf_object = new TIGR::Foundation;
      $tf_object->TIGR_GetOptions();
      my @error_arr = ();
      my $error_arr_ref = \@error_arr;
      # create new ConfigFile object
      my $cf = new TIGR::ConfigFile("test_conf", $tf_object,
                  $error_arr_ref);

      if(defined $cf) {
         # get the sections in the config file
         my $section_ref = $cf->getSections();
         my @sections = @$section_ref;
         my $section = undef;
         my @options = ();
         my $option = undef;
         my $option_val = undef;
         my $options_ref = {};

         foreach $section (@sections) {
            # get the options for each section
            $options_ref = $cf->getOptionNames($section);
            @options = @$options_ref;

            foreach $option (@options) {
               # get the value for the option
               $option_val = $cf->getOption($option, $section, "default");
               print("the value of option $option for section $section is ".
                  "$option_val\n");
            }
         }
      }
   }

=cut

}
1;

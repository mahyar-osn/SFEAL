#!/usr/bin/perl

$cmiss_root=$ENV{CMISS_EXECUTABLE};  # define cmiss root
$cmgui_root=$ENV{CMGUI_2_6_2};  # define cmgui root
$ip2ex_cm=$ARGV[0];
$cmgui_file=$ARGV[1];
$show_cmgui=$ARGV[2];

$comfile = $ip2ex_cm;
$cmguifile = $cmgui_file;

system "$cmiss_root $comfile";
if($show_cmgui =~ 'show') {system "$cmgui_root $cmguifile";}


# PSD statistics processor, written by wjm34 24/Jun/2011

proc psl-stats {inpfile} {
	#Read psl file
	set data  ""
	set dcoldata ""
	set n 0
	set pi 3.1415926535
	set fname [open $inpfile r]
	set pri_index 13
	while { [gets $fname line] != -1 } {
		set fields [split $line ","]
		if {$n == 0} {
		   # for {set i 0} {$i < [llength $line]} {incr i} {
		   #     if {[regexp "Primary Diameter" [lindex $line $i]]} {
           #     set pri_index $i
		   #     }
		   # }
		}  else {
			lappend data [lindex $fields $pri_index]
		}
		incr n
    }
	close $fname
	
    set dpri_mean [calc $data]
    return $dpri_mean
}	
	
#Subroutine to calculate statistics
proc calc {inp} {
  set l [llength $inp]
  
  #Calculate arithmetic and geometric means
  set sum 0
  set sumsq 0
  set prod 1
  set p [expr 1.0/$l]
  for {set i 0} {$i < $l} {incr i} {
    set sum [expr $sum + [lindex $inp $i]]
    set sumsq [expr $sumsq + [lindex $inp $i] * [lindex $inp $i]]
    set prod [expr $prod * pow([lindex $inp $i],$p)]
  }
  set amean [expr $sum / $l]
  set a2mean [expr $sumsq / $l]
  set gmean $prod

  return [list $amean $a2mean $l]
} 

# Calculate data from files
set test_finite [psl-stats silane-finite-psl(0.08s).csv]
set test_finite_mean [lindex $test_finite 0]
set test_finite_conf [expr sqrt(([lindex $test_finite 1] \
                                 - $test_finite_mean * $test_finite_mean) \
                                / [lindex $test_finite 2]) * 3.29]
set test_spherical [psl-stats silane-spherical-psl(0.08s).csv]
set test_spherical_mean [lindex $test_spherical 0]
set test_spherical_conf [expr sqrt(([lindex $test_spherical 1] \
                                   - $test_spherical_mean * $test_spherical_mean) \
                                  / [lindex $test_spherical 2]) * 3.29]
set test_nosinter [psl-stats silane-nosinter-psl(0.08s).csv]
set test_nosinter_mean [lindex $test_nosinter 0]
set test_nosinter_conf [expr sqrt(([lindex $test_nosinter 1] \
                                  - $test_nosinter_mean * $test_nosinter_mean) \
                                 / [lindex $test_nosinter 2]) * 3.29]

# Test if it is working!

set val_finite 26.3
set val_spherical 42.6
set val_nosinter 0.49

if { abs($test_finite_mean - $val_finite) < 0.3 } {
    puts "Particle diameter $test_finite_mean. Passed! (conf: $test_finite_conf)"
} else {
    puts "TEST FAILED: diameter was $test_finite_mean, $val_finite expected.(conf: $test_finite_conf)"
    exit 1
}

if { abs($test_nosinter_mean - $val_nosinter) < 0.01 } {
    puts "Particle diameter $test_nosinter_mean. Passed! (conf: $test_nosinter_conf)" 
} else {
    puts "TEST FAILED: diameter was $test_nosinter_mean, $val_nosinter expected. (conf: $test_nosinter_conf)" 
    exit 1
}


if { abs($test_spherical_mean - $val_spherical) < 1.5 } {
    puts "Particle diameter $test_spherical_mean. Passed! (conf: $test_spherical_conf)"
} else {
    puts "TEST FAILED: diameter was $test_spherical_mean, $val_spherical expected. (conf: $test_spherical_conf)"
    exit 1
}

puts "All tests passed."
exit 0

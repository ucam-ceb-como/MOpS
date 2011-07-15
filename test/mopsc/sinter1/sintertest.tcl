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
  set prod 1
  set p [expr 1.0/$l]
  for {set i 0} {$i < $l} {incr i} {
    set sum [expr $sum + [lindex $inp $i]]
	set prod [expr $prod * pow([lindex $inp $i],$p)]
	}
  set amean [expr $sum / $l]
  set gmean $prod

  return $amean
} 

# Calculate data from files
set test_finite [psl-stats silane-finite-psl(0.08s).csv]
set test_spherical [psl-stats silane-spherical-psl(0.08s).csv]
set test_nosinter [psl-stats silane-nosinter-psl(0.08s).csv]

# Test if it is working!
# allow 3% tolerance
set utol 1.03
set ltol 0.97

set val_finite 26.3
set val_spherical 42.6
set val_nosinter 0.49

if {$test_finite < $utol*$val_finite && $test_finite > $ltol*$val_finite} {
    puts "Particle diameter $test_finite. Passed!"
} else {
    puts "TEST FAILED: diameter was $test_finite, $val_finite expected."
    exit 1
}

if {$test_nosinter < $utol*$val_nosinter && $test_nosinter > $ltol*$val_nosinter} {
    puts "Particle diameter $test_nosinter. Passed!"
} else {
    puts "TEST FAILED: diameter was $test_nosinter, $val_nosinter expected."
    exit 1
}


if {$test_spherical < $utol*$val_spherical && $test_spherical > $ltol*$val_spherical} {
    puts "Particle diameter $test_spherical. Passed!"
} else {
    puts "TEST FAILED: diameter was $test_spherical, $val_spherical expected."
    exit 1
}

puts "All tests passed."
exit 0

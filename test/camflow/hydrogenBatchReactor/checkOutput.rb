#!/usr/bin/ruby -w
# Check H2 Batch Reactor Output.
# Laurence R. McGlashan

old = Array.new
new = Array.new;

# Open the newly simulated solution file.
File.open("profile.dat", "r") do |infile|
    while (line = infile.gets)
        new.clear
        new += line.split
    end
end

# Open the file containing the reference solution.
File.open("profileOriginal.dat", "r") do |infile|
    while (line = infile.gets)
        old.clear
        old += line.split
    end
end

# Compare the strings of fields 3 to 10. They should match.
(3..10).each do |i|
    if (old[i] != new[i])
      puts "Old value was #{old[i]}, but new value is #{new[i]}"
      puts "**************************"
      puts "****** TEST FAILURE ******"
      puts "**************************"
      return 1
    end
end

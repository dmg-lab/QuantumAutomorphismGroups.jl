import FileWatching: watch_file 

In_loc =  try ARGS[1] catch _ "input.txt" end

watch_file(In_loc);
println("Something changed in $In_loc")


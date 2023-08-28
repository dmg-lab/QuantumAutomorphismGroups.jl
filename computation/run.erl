
-module(run).

-export([
         main/0,
        read/0,
        decode/1,
        loop/1,
        start/0,
        init_worker/0,
        watch_file/1,
        compute/1,
        stop/1
        ]).


main() ->
     %file:set_cwd("../computation/"),
    read(), 
    %Cmd = "julia ./run.jl",
    %Path = "uniform_matroid/r3n4",
    %Code = "uniform_matroid(3,4)",
    {ok, Cwd} = file:get_cwd(),
    io:format("Current working directory: ~s~n", [Cwd]),
    io:format("Main Pid: ~p~n", [self()]),
    start_file_watcher("input.txt").

    %spawn(fun() -> start_monitor() end),
    %Worker_Pid = start(),
    %compute(Worker_Pid).



%start_monitor() ->
%    register(monitor, self()),
%    Cmd = "mkdir -p data",
%    os:cmd(Cmd),
%    io:format("Monitoring...~n"),
%    monitor_loop().
%
%
%monitor_loop() ->
%    receive
%        {worker, Msg} ->
%            io:format("Worker finished successfully: ~p~n~s~n", [worker, Msg]),
%            monitor_loop();
%        _ ->
%            io:format("Something went wrong~n"),
%            monitor_loop()
%    end.

start_file_watcher(File) ->
    Watcher_Pid = spawn(?MODULE, watch_file, [File]),
    Watcher_Pid.
    
watch_file(File) ->
    {Julia_Pid, Monitor} = spawn_monitor(fun() -> os:cmd("julia ./utils/filewatcher.jl " ++ File) end),
    monitor(process, Julia_Pid),
    receive
        {'DOWN', Monitor, process, _ , Info} ->
            io:format("File changed: ~p~n", [Info]),
            watch_file(File)
    end.





read() ->
    %{ok, [Input]} = io:fread("Enter the name of the file: ", "~s"),
    Input = "input.txt",
    Otp = os:cmd("julia ./utils/read.jl " ++ Input),
    Arr = string:tokens(Otp, "\n"),
    Sep = fun (Line) -> string:tokens(Line, ":") end,
    LineList = [Sep(L) || L <- Arr],
    io:format("~p~n", [LineList]).
    %io:format("~s~n", [FirstLine]).

decode(Data) ->
    {_, Msg} = Data,
    Msg.


loop(Port) ->
    receive
	{call, Caller, Msg} ->
	    Port ! {self(), {command, Msg}},
	    receive
		{Port, {data, Data}} ->
		    Caller ! {self(), decode(Data)}
	    end,
	    loop(Port);
	stop ->
	    Port ! {self(), close},
	    receive
		{Port, closed} ->
		    exit(normal)
	    end;
	{'EXIT', Port, _} ->
	    exit(port_terminated)
    end.


start() ->
    Worker_Pid = spawn(?MODULE, init_worker, []),
    io:format("Created process: ~p~n", [Worker_Pid]),
    Worker_Pid.

init_worker() ->
    Port = open_port({spawn, "julia"}, [{line, 256}, use_stdio,  eof, exit_status]),
    loop(Port).

compute(Worker_Pid) ->
    Worker_Pid ! {call, self(), "x=2; x+4 \n"},
    receive
	    {Worker_Pid, Result} ->
            io:format("Result: ~p~n", [Result])
    end.

stop(Worker_Pid) ->
    Worker_Pid ! stop.


read_tst := proc(tst_number::posint, $)
	local str;

	str := cat("tst/tst", 
			    convert(tst_number, string), ".mm"); 
	try 	
		read(str);
		print("Prime characteristic:", m);
		print("f: ", f);
		print("t: ", t);
		print("b: ", b);
	catch "unable to read `%1`":
		error "test number 1% does not exist", tst_number;
	end try;

end proc;


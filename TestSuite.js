 function TestSuite(about){
	var passed_tests = 0;
	var run_tests = 0;
	var about = (about===undefined?"Empty info about the test":about.toString());
	var test_info = "";
	this.run_test = function(computed_value, expected_value, info){
		run_tests++;
		if (info==undefined) {info="";};
		var error_message = "";
		error_message += "Test #" + run_tests + " failed: " + info +
                        " \nComputed: " + computed_value +
                        " \nExpected: " + expected_value;
		if (computed_value==expected_value) {passed_tests++;} 
		else {alert(error_message);};
	};
	var to_String = function(){
		test_info="";
		test_info += "Results: Passed " + passed_tests + " of " + run_tests + ".";
	};
	this.alert_results = function(){
		to_String();
		var message = about + "\n" + test_info;
		alert(message);
	}
	this.results_on_console = function(){
		to_String();
		var message = about + "\n" + test_info;
		console.log(message);
	}
	this.results_on_document = function(){
		to_String();
		var message = about + "<BR>" + test_info;
		document.write(message);
	}
}

 function TestSuite(about){
	/*	This object helps to test different parts of the program.
	 *	Use:
	 *	    - Create a TestSuite object:    ts = new TestSuite("Info about this group of tests");
     *      - Run a test:                   ts.run_test(Computed value, Expected value, "error message");
     *      - Check results:                ts.results_on_console();
	 */

    //variables
	var passed_tests = 0;
	var run_tests = 0;
	var about = (about===undefined?"Empty info about the test":about.toString());
	var test_info = "";


    //methods
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

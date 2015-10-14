/*  Matrix library
 *	Author: Jaime Pérez Aparicio
 *	email: 19.jaime.91@gmail.com
 *	license: GPL
 *
 *
 */



function Matrix(rows, cols, data){
	// On creatiion:
    this._rows=(rows==undefined?1:rows);
    this._cols=(cols==undefined?1:cols);
	this.create=function(){
        __data=new Array(this._rows);
        for(var row=0; row<this._rows; row++){
            __data[row]=new Array(this._cols);
        }
		return __data;
	}
	this.data_correct=function(){
		if (data===undefined){
			return false;
		} 
		return true;
	}
	this.prepare=function(_data){
		var output = [];
		for (var row=0; row<_data.length; row++){
			output.push(data[row].slice(0));	
		}
		return output;
	}
    this._data=(this.data_correct()?this.prepare(data):this.create());
    this.show=function(){
        document.write(this._data);
    }
    this.getRows=function(){
        return this._rows;
    }
    this.getCols=function(){
        return this._cols;
    }
    this.get=function(row, col){
		try{
			output = this._data[row][col];
		}catch(exception){
			output = undefined;
			alert("Bad matrix!");}
		return output;
    }
    this.set=function(row, col, value){
        this._data[row][col] = value;   
    }
    this.set_zeros=function(){
        for(var row=0; row<this.getRows(); row++){
            for(var col=0; col<this.getCols(); col++){
                this._data[row][col] = 0;
            }
        }
    }
    this.set_homogeneus_value=function(value){
        for(var row=0; row<this.getRows(); row++){
            for(var col=0; col<this.getCols(); col++){
                this._data[row][col] = value;
            }
        }
    }
    this.minor=function(_row,_col){
        matrix = new Matrix(this.getRows()-1,this.getCols()-1);
        var __row=_row;
        var __col=_col;
        
            for(var row=0; row<this.getRows(); row++){
                if(row<_row){
                    __row=row;
                }else if(row>_row){
                    __row=row-1;
                };
                for(var col=0; col<this.getCols(); col++){
                    if(col<_col){
                        __col=col;
                        matrix.set(__row,__col,this.get(row,col));
                    }else if(col>_col){
                        __col=col-1;
                        matrix.set(__row,__col,this.get(row,col));
                    };
                }
        }
        return matrix;
    }	
    this.determinant=function(){
		/*
		 *
		 *	Computes determinant with a recursive method.
		 *
		 */
        var total = 0;
        if(this.getRows() != this.getCols()) return undefined;
		if(this.getRows() == 2){
			var result = this.get(0,0)*this.get(1,1);
			result -= this.get(0,1)*this.get(1,0);
			return result;
		}
		for (var row=0; row<this.getRows(); row++){
			var minor = this.minor(row,0);
			total += minor.determinant()*(-1)^row;
		}
        return total;
    }
    this.trace=function(){
        if (this.getRows()!=this.getCols()){return undefined;};
        var trace=0;
        for(var index=0; index<this.getCols(); index++){
            trace+=this.get(index,index);
        }
        return trace;
    }
	this.equals=function(matrix){
		if (this.id==matrix.id) return true;
		condition = false;
		if(typeof(matrix) == Matrix){
			if (this.get_rows()==matrix.get_rows()){
				if (this.get_cols()==matrix.get_cols()){
					condition = true;
				}
			}
		}
		if (!condition) return false;
		for (var row=0; row<this.getRows(); row++){
			for (var col=0; col<this.getCols(); col++){
				if (this.get(row,col) != matrix.get(row, col)) return false;
			}
		}
		return condition;
	}
	this.clone=function(){
		var rows = this.getRows();
		var cols = this.getCols();
		var data = this._data;
		var output = new Matrix(rows, cols, data);
		return output;
	}
	this.traspose=function(){
		return;
	}
}

function Tests(){
    ts = new TestSuite("Matrix class test");
    mat = new Matrix(2,2);
    mat.set_zeros();
    ts.run_test(mat.get(0,0),0,"Set zeros method failed!");
    mat.set_homogeneus_value(1  );
    ts.run_test(mat.get(0,0),1,"Set homogeneus method failed!");  
    mat2 = mat.minor(1,1);
	mat3=new Matrix(2,2);
	mat3.set_zeros();
	ts.run_test(mat.equals(mat3),true,"equals method failed!");
	mat4 = new Matrix();
	mat4.set_zeros();
	ts.run_test(mat2.equals(mat4),true,"minor method failed!");
	ts.run_test(mat.determinant(),0,"Determinant computation failed!");
	ts.run_test(mat.trace(),2,"Trace computation method failed!");
	mat3 = mat.clone();
	ts.run_test(mat.equals(mat3),true,"clone method failed!");
	v1 = [1,2,2];
	v2 = v1.slice(1);
	ts.run_test(v1[1]==v2[0], true,"Slice doesn't work as expected!");
	v2[0] = 5;
	ts.run_test(v1[1]==v2[0], false,"Slice only references!");
	mat3.set(0,0,5);
	ts.run_test(mat.equals(mat3),false,"clone method failed! Both matrix changed! Check creation method (with input data).");
	ts.results_on_console();
}

new Tests();

module LinearCodes :
	sig
		val print_matrix : int array -> unit
		val getCode : int array -> int -> int
		val linear_even : int -> int
		val systematic_matrix : int -> int -> int array -> int array
		val matrix_transp : int array -> int array
		val systematic_ctrl_matrix : int -> int -> int array -> int array
		val getSyndrome : int array -> int -> int
		val compute_code : int array -> int -> int
		val compute_correction : int array -> int -> int
	end =
	struct

	(* Display a bit-matrix *)
	let print_matrix m =
		(* Bits to string *)
		let str_m = Array.map BitsTools.bin_of_int m in
		(* Find the longest string length - keep string sizes in memory*)
		let n = Array.length str_m in
		let str_s = Array.make n 0 in
		let max_l = ref 0 in
		for i = 0 to (n-1) do
			let cur_l = (String.length str_m.(i)) in
			str_s.(i)<-cur_l;
			if (cur_l > !max_l) then max_l := cur_l;
		done;
		(* Make lists with zeroes *)
		let rec makeZeroString x =
			if x = 0 then ""
			else (String.concat "" ["0";(makeZeroString (x-1))]) in
		(* Pad everyone - tidious *)
		for i = 0 to (n-1) do
			str_m.(i)<-(String.concat "" [(makeZeroString (!max_l-str_s.(i)));str_m.(i)]);
		done;
		(* Display the matrix *)
		for i = 0 to (!max_l-1) do
			for j = 0 to (n-1) do
				print_char str_m.(j).[i];
				print_string " ";
			done;
			print_string "\n";
		done;;

	(* Naive bit-matrix implementation *)
	let getCode matrix vect =
		(* Multiply and sum the vectors *)
		let multiply a b =
			let cur_a = ref a in
			let cur_b = ref b in
			let i = ref 0 in
			let ret = ref 0 in
			while (!cur_a <> 0 && !cur_b <> 0) do
				ret := !ret + ((!cur_a land 1) * (!cur_b land 1));
				(* Switch to next bit *)
				cur_a := (!cur_a lsr 1);
				cur_b := (!cur_b lsr 1);
				i := !i + 1
			done;(!ret mod 2) in
		(* Actual multiplication *)
		let n = Array.length matrix in
		let ret = ref 0 in
		for i = 0 to (n-1) do
			ret := !ret + ((multiply vect matrix.(i)) lsl (n-i-1))
		done;!ret;;

	(* Linearized even code *)
	let linear_even a =
		(* Get size of int *)
		let size = (int_of_float (BitsTools.sizeof a)) in
		(* Make the matrix *)
		let matrix = Array.make (size+2) 0 in
		matrix.(0)<-(1 lsl size);
		for i = 1 to (size) do
			matrix.(i)<-matrix.(i-1)/2;
		done;
		matrix.(size+1)<-((matrix.(0)*2)-1);
		(* Apply it on code *)
		getCode matrix a;;

	(* Create systematic code matrix *)
	let systematic_matrix n k b =
		(* Make identity matrix *)
		let matrix = Array.make n 0 in
		matrix.(0)<-(1 lsl (k-1));
		for i = 1 to (k-1) do
			matrix.(i)<-matrix.(i-1)/2;
		done;
		(* Complete with B matrix *)
		let rest = Array.length b in
		for i = 0 to (rest-1) do
			matrix.(k+i)<-b.(i);
		done;matrix;;

	(* Bit-matrix transposition *)
	let matrix_transp matrix =
		let n = Array.length matrix in
		(* Get the max size - Unefficient *)
		let size = ref (BitsTools.sizeof matrix.(0)) in
		for i = 1 to (n-1) do
			let cur_size = (BitsTools.sizeof matrix.(i)) in
			if cur_size > !size then size := cur_size;
		done;
		(* Do the transposition *)
		let size_i = (int_of_float !size) in
		let ret = Array.make size_i 0 in
		for i = 0 to (n-1) do
			for j = 0 to (size_i-1) do
				ret.(j) <- (((matrix.(i) lsr (size_i-j-1)) land 1) lsl (n-i-1)) + ret.(j);
			done;
		done;ret;;
			
	(* Create systematic control matrix *)
	let systematic_ctrl_matrix n k b =
		let ret = Array.make n 0 in
		(* Do the left side of the matrix *)
		let b_t = matrix_transp b in
		let b_s = Array.length b in
		for i = 0 to (b_s-1) do
			ret.(i) <- b_t.(i);
		done;
		(* Do the idenity *)
		ret.(n-k)<-(1 lsl (k-1));
		for i = (n-k+1) to (n-1) do
			ret.(i)<-(ret.(i-1) lsr 2);
		done;ret;;

	(* Get syndrome *)
	let getSyndrome h a = getCode (matrix_transp h) a;;

	(* Calculate code of a bit-sequence *)
	let compute_code g a = getCode g a;;

	(* Find error in a bit sequence and correct it *)
	let compute_correction h a =
		(* Get the syndrome *)
		let s = getSyndrome h a in
		(* Find the error location *)
		let error = ref (-1) in
		for i = 0 to ((Array.length h)-1) do
			if (h.(i) = s) then error := i;
		done;
		(* Fix the error *)
		if (!error = (-1)) then a
		else
			(let size = BitsTools.sizeof a in
			a lxor (1 lsl ((int_of_float size) - !error)));;

	end

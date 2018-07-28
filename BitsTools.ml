(* Libraries *)
#load "unix.cma";;
open Unix;;

module BitsTools :
	sig
		val or_gate : int -> int -> int
		val and_gate : int -> int -> int
		val xor_gate : int -> int -> int
		val not_gate : int -> int
		val right_shift : int -> int -> int
		val left_shift : int -> int -> int
		val bin_of_int : int -> string
		val pow : int -> int -> int
		val hamming_weight_naive : int -> int
		val hamming_weight_sparse : int -> int
		val lg : float -> float
		val sizeof : int -> float
		val alignment : float -> int ref -> int
		val hamming_weight_dense : int -> int
		val hamming_weight_parallel : int -> int
		val hamming_weight_wp3 : int -> int
		val time : ('a -> 'b) -> 'a -> unit
		val benchmark : int -> int -> unit
		val hamming_distance : int -> int -> int
		val hamming_weight : int -> int
	end =
	struct

	(* Structure *)
	type bit_array = int;;

	(* OR Gate *)
	let or_gate a b =  a lor b;;

	(* AND Gate *)
	let and_gate a b = a land b;;

	(* XOR Gate *)
	let xor_gate a b = a lxor b;;

	(* NOT Gate *)
	let not_gate a = lnot a;;

	(* Shuffle *)
	let right_shift a dir = a lsr dir;;
	let left_shift a dir = a lsl dir;;

	(* Taken from https://rosettacode.org/wiki/Binary_digits#OCaml *)
	let bin_of_int d =
		if d < 0 then invalid_arg "bin_of_int" else
		if d = 0 then "0" else
		let rec aux acc d =
		  if d = 0 then acc else
		  aux (string_of_int (d land 1) :: acc) (d lsr 1)
		in
		String.concat "" (aux [] d);;

	(* Integer efficient power funtion *)
	let rec pow a = function
		| 0 -> 1
		| 1 -> a
		| n ->
			let b = pow a (n / 2) in
			b * b * (if n mod 2 = 0 then 1 else a);;

	(* Hamming Weight *)
	(* Iterate everything *)
	let rec hamming_weight_naive a =
		if a = 0 then 0
		else (a land 1) + (hamming_weight_naive (a lsr 1));;
	(* Iterate as much time as there are 1-bits *)
	let rec hamming_weight_sparse a =
		if a = 0 then 0
		else 1 + (hamming_weight_sparse (a land (a-1)));;
	(* Get the size of the int *)
	let lg x = ((log x) /. (log 2.));;
	let sizeof x = floor (lg (float_of_int x)) +. 1.;;
	let alignment size x =
		(* Get the closest power of two for size *)
		let next = (2.**(ceil (lg size))) in
		(* Align the number *)
		x := (!x lsl (int_of_float (next -. size)));
		(* Return the aligned basis *)
		(int_of_float (2.**next));;
	(* Dense-ones: Only iterates as many times as there are zero-bits in the integer.*)
	let hamming_weight_dense a =
		let rec hamming_weight_dense_aux count a =
			if a = 0 then count
			else (hamming_weight_dense_aux (count-1) (a land (a-1))) in
		(hamming_weight_dense_aux (8*8-1) (lnot a));;
	(* Nifty parallel counting. *)
	let rec hamming_weight_parallel a =
		let n = ref a in
		(* Split the number if it's too large (more than 32 bits long) *)
		let size = sizeof a in
		if (size > 32.) then
			let cut = (floor (size /. 2.)) in
			let a_low = (a land ((int_of_float (2.**cut))-1)) in
			let a_high = (a lsr (int_of_float (size -. cut))) in
			((hamming_weight_parallel a_low) + (hamming_weight_parallel a_high))
		else
		(* Get re-alignment goodies *)
		let align =  alignment size n in
		(* Do the magic *)
		let m1 = align / 3 in (* Binary 0101010101... *)
		let m2 = align / 5 in (* Binary 0011001100... *)
		let m4 = align / 17 in (* Binary 000011110000 ... *)
		n := (!n land m1) + ((!n lsr 1) land m1);
		n := (!n land m2) + ((!n lsr 2) land m2);
		n := (!n land m4) + ((!n lsr 4) land m4);
		!n mod 255;;
	(* WP3 - Nifty Revised *)
	let rec hamming_weight_wp3 a =
		let n = ref a in
		(* Split the number if it's too large (more than 32 bits long) *)
		let size = sizeof a in
		if (size > 32.) then
			let cut = (floor (size /. 2.)) in
			let a_low = (a land ((int_of_float (2.**cut))-1)) in
			let a_high = (a lsr (int_of_float (size -. cut))) in
			((hamming_weight_parallel a_low) + (hamming_weight_parallel a_high))
		else
		(* Get re-alignment goodies *)
		let align =  alignment size n in
		(* Do the magic *)
		let m1 = align / 3 in (* Binary 0101010101... *)
		let m2 = align / 5 in (* Binary 0011001100... *)
		let m4 = align / 17 in (* Binary 000011110000 ... *)
		n := !n - ((!n lsr 1) land m1); (* Put count of each 2 bits into those 2 bits *)
		n := (!n land m2) + ((!n lsr 2) land m2); (* Put count of each 4 bits into those 4 bits *)
		n := (!n + (!n lsr 4)) land m4; (* Put count of each 8 bits into those 8 bits *)
		!n mod 255;;

		(* Execution time of a function *)
		let time f x =
			let start = Unix.gettimeofday ()
			in let res = f x
			in let stop = Unix.gettimeofday ()
			in let () = Printf.printf "%fs%!" (stop -. start)
			in res; ();;

		(* Benchmark Hamming Weight Algorithms *)
		let benchmark bits_length iter =
			(* Get int limit *)
			let limit = (int_of_float (2.**(float_of_int bits_length))) in
			for i = 0 to iter do
				(* Get random int *)
				let a = (Random.int (limit+1)) in
				(* Get size *)
				print_string "Size = ";
				print_float (sizeof a);
				(* Naive algorithm *)
				print_string " | Naive = ";
				time hamming_weight_naive a;
				(* Sparse algorithm *)
				print_string " | Sparse = ";
				time hamming_weight_sparse a;
				(* Dense algorithm *)
				print_string " | Dense = ";
				time hamming_weight_dense a;
				(* Parallel algorithm *)
				print_string " | Parallel = ";
				time hamming_weight_parallel a;
				(* WP3 algorithm *)
				print_string " | WP3 = ";
				time hamming_weight_wp3 a;
				(* Skip line *)
				print_string "\n";
			done;;

		(* Hamming Weight *)
		let hamming_weight a = hamming_weight_naive a;;

		(* Hamming distance *)
		let hamming_distance a b = hamming_weight_naive (b - a);;

	end

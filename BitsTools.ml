(* Libraries *)
#load "unix.cma";;
open Unix;;

module BitsTools :
	sig
		type bit_array = int
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
		val hamming_weight_dense : int -> int
		val hamming_weight_parallel : int -> int -> int
		val hamming_weight_wp3 : int -> int -> int
		val time : ('a -> 'b) -> 'a -> unit
		val benchmark : int -> unit
		val hamming_weight : int -> int -> int
		val hamming_distance : int -> int -> int -> int
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
	(* Dense-ones: Only iterates as many times as there are zero-bits in the integer.*)
	let hamming_weight_dense a =
		let rec hamming_weight_dense_aux count a =
			if a = 0 then count
			else (hamming_weight_dense_aux (count-1) (a land (a-1))) in
		(hamming_weight_dense_aux (8*8-1) (lnot a));;
	(* Nifty parallel counting. *)
	let rec hamming_weight_parallel a size =
		let n = ref a in
		(* Do the magic *)
		let size_int = (1 lsl size) in
		let m1 = size_int / 3 in (* Binary 0101010101... *)
		let m2 = size_int / 5 in (* Binary 0011001100... *)
		let m4 = size_int / 17 in (* Binary 000011110000 ... *)
		n := (!n land m1) + ((!n lsr 1) land m1);
		n := (!n land m2) + ((!n lsr 2) land m2);
		n := (!n land m4) + ((!n lsr 4) land m4);
		!n mod 255;;
	(* WP3 - Nifty Revised *)
	let rec hamming_weight_wp3 a size =
		let n = ref a in
		(* Do the magic *)
		let size_int = (1 lsl size) in
		let m1 = size_int / 3 in (* Binary 0101010101... *)
		let m2 = size_int / 5 in (* Binary 0011001100... *)
		let m4 = size_int / 17 in (* Binary 000011110000 ... *)
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
		let benchmark iter =
			(* Get int limit *)
			for i = 0 to iter do
				(* Get random int *)
				let a = (Random.bits ()) in
				if (a <= 4294967296) then (* Only do 32 bits *)
					(* Naive algorithm *)
					(print_string " | Naive = ";
					time hamming_weight_naive a;
					(* Sparse algorithm *)
					print_string " | Sparse = ";
					time hamming_weight_sparse a;
					(* Dense algorithm *)
					print_string " | Dense = ";
					time hamming_weight_dense a;
					(* Parallel algorithm *)
					print_string " | Parallel = ";
					time (hamming_weight_parallel a) 32;
					(* WP3 algorithm *)
					print_string " | WP3 = ";
					time (hamming_weight_wp3 a) 32;
					(* Skip line *)
					print_string "\n";);
			done;;

		(* Hamming Weight *)
		let hamming_weight a size = hamming_weight_parallel a size;;

		(* Hamming distance *)
		let hamming_distance a b size = hamming_weight_parallel (b - a) size;;
	end

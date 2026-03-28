all:
	cargo run 

gitaddall:
	git add src crates examples tests

loc:
	find src tests crates examples -name '*.rs' | xargs wc -l

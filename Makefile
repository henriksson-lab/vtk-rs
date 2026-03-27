all:
	cargo run 

gitaddall:
	git add src crates

loc:
	find src tests crates -name '*.rs' | xargs wc -l

all:
	cargo run 

gitaddall:
	git add src ui assets crates

loc:
	find src tests crates -name '*.rs' | xargs wc -l

clean:cleancrap cleanexe

cleancrap:
	@echo "Cleaning crap..."
	@find . -name "*~" -exec rm -rf {} \;
	@find . -name "#" -exec rm -rf {} \;

cleanexe:
	@echo "Cleaning byte code..."
	@rm *.pyc

commit:
	@echo "Commiting..."
	@-git commit -am "Commit"
	@-git push origin master

pull:
	@-git reset --hard HEAD
	@-git pull


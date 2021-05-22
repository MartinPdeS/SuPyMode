
upload:
	python3.8 -m twine upload --password ${PipPassword} --username ${PipUsername} --repository pypi ./dist/*

read:
	google-chrome docs/build/html/index.html

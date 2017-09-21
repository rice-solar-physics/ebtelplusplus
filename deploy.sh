scons
cd docs && make && cd ..
# ghp-import -m "Build docs" -b gh-pages docs/html
# git push -fq https://${GH_TOKEN}@github.com/${TRAVIS_REPO_SLUG} gh-pages
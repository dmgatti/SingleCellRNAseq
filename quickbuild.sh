# just build the jekyll website locally from markdown already generated
bundle config set --local path .vendor/bundle && bundle install && bundle update && bundle exec jekyll build --incremental

echo "On a mac you can show the website by executing the command:"
echo "open _site/index.html"

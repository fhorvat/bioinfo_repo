### TUTORIAL: https://product.hubspot.com/blog/git-and-github-tutorial-for-beginners
# set message
MESSAGE=""

### go to dir, add all new changes, commit
# go to main dir with scripts
cd /common/WORK/fhorvat/Projekti/Svoboda/scripts

# add files to staging environment (--all will take in account removed dirs/files)
git add --all .

# create a commit
git commit -m "$MESSAGE" -m "`$commitDate`"


#### create a new branch
## optionally create a new branch
#git checkout -b <my branch name>
#
## switch back to master branch
##git checkout master
#
## delete branch
##git branch -d <my branch name>


### push
# push a master branch to GitHub (username: fhorvat)
git push bioinfo master

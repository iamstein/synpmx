## Vibe coding notes

This is a sort of journal of my (Andy's) vibe coding experience as some things to be careful of.  

It's tempting to generate code quickly without truly understanding all that is happening.  you get working code right away!  but this also caused some issue.  Here were some.

- I didn't initially realize that AVATAR, when there was only one patient at a dose level receiving a dose at a certain time, the algortihm wouldn't blend patients, it would just add noise to teh data.  Similarly, if there was only one patient 
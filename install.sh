#!/bin/bash
path=`pwd`
mkdir logs

echo "Checking if Bowtie is in your \$PATH...";
if command -v bowtie 1> /dev/null; then
    echo "Bowtie exists.";
	bowtie --version
else
    echo "Bowtie was not found in your \$PATH. I will try to download it for you from the GIT repository.";
    
    git clone https://github.com/BenLangmead/bowtie.git
    cd bowtie
	echo -n "Building Bowtie"
    make -s 2>$path/logs/bowtie_build.log
	echo ", done."

    if [ ! -d ~/bin ];then
        echo "Creating bin folder in the home directory.";
        mkdir ~/bin
        echo "The local bin folder will now be included in your \$PATH (see $HOME/.profile)."

		if [ -e ~/.profile ];then
			echo "\n\n# This line has been inserted by the QUASI Pipeline\nexport PATH=$PATH:$HOME/bin" >> ~/.profile;
			echo "$(tput setaf 1;tput bold)The path has been added to $HOME/.profile. Please execute the following after this script has ended: source ~/.profile$(tput sgr0)"

		elif [ -e ~/.bash_profile ];then
			echo "\n\n# This line has been inserted by the QUASI Pipeline\nexport PATH=$PATH:$HOME/bin" >> ~/.bash_profile
			echo "$(tput setaf 1;tput bold)The path has been added to $HOME/.bash_profile. Please execute the following after this script has ended: source ~/.bash_profile$(tput sgr0)"

		else
			echo "\n\n# This line has been inserted by the QUASI Pipeline\nexport PATH=$PATH:$HOME/bin" >> ~/.profile
			echo "$(tput setaf 1;tput bold)The path has been added to $HOME/.profile. Please execute the following after this script has ended: source ~/.profile$(tput sgr0)"
		fi

    else
        echo "Found bin in the home directory. Installing Bowtie there.";
    fi
	echo "Moving the bowtie executables to $HOME/bin.";
	find . -perm /a+x -type f -name "bowtie*" -exec cp {} ~/bin \;
	cd ..
fi

if [ ! -e qa ];then
    echo -n "I will now create the quality assessment tool"
    make -s qa 2>$path/logs/qa_build.log
	echo ", done."
fi

if [ ! -e count ];then
    echo -n "I will now create the quantification tool"
    make -s count 2>$path/logs/count_build.log
	echo ", done."
fi

echo "Creating the reference folder in your home directory. Place your reference fasta files and Bowtie/BWA index files there."
mkdir ~/references

echo "Finished installing. You may now run the pipeline.";

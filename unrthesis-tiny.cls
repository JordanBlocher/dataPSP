%%
%% File:   unrthesis.cls
%%

% Building off of report class
\ProvidesClass{unrthesis}
\LoadClass{report}

% New commands
\newcommand{\@titlesize}{\fontsize{18}{24}\selectfont}
\newcommand{\@subtitlesize}{\fontsize{12}{24}\selectfont}
\newcommand{\@fontsize}{\fontsize{12}{24}}
\newcommand{\@authorsize}{\normalsize\doublespacing} 
\newcommand{\@advisor}{Dr. Eric Olson}
\newcommand{\@masterstype}{Science}
\newcommand{\@degreefield}{Applied Math}

\newcommand{\@conferralmonth}{May}
\newcommand{\@conferralyear}{2016}


% Switch between Masters and PhD
\newif\ifphd
\DeclareOption{phd}{\phdtrue}
\DeclareOption{master}{\phdfalse}
\ExecuteOptions{master}


% Create list of committee members
\newcommand\addcommitteemember[3]{%
  \write\@auxout{\noexpand\@writefile{members}{#1, #2 -- #3\noexpand\linebreak ~\noexpand\linebreak}}}
\newcommand\addgraduatedean[2]{%
  \write\@auxout{\noexpand\@writefile{members}{#1, #2 -- Dean, Graduate School}}}
\newcommand\listofmembers{%
    \@starttoc{members}%
}
\addcommitteemember{Eric Olson}{Ph.\ D.}{Advisor}                           %
\addcommitteemember{Bruce Blackadar}{Ph.\ D.}{Committee Member}                %
\addcommitteemember{Paul J. Hurtado}{Ph.\ D.}{Committee Member} % the committee
\addcommitteemember{Fred Harris}{Ph.\ D.}{Graduate School Representative}
\addgraduatedean{David Zeh}{Ph.\ D.}   

\newcommand{\realsinglespacing}{%
  \let\CS=\@currsize\renewcommand{\baselinestretch}{1.0}\CS}
\newcommand{\singlespacing}{%
  \let\CS=\@currsize\renewcommand{\baselinestretch}{1.208333}\CS}
\newcommand{\singlespacingplus}{%
  \let\CS=\@currsize\renewcommand{\baselinestretch}{1.5}\CS}
\newcommand{\doublespacing}{%
  \let\CS=\@currsize\renewcommand{\baselinestretch}{2}\CS}
\newcommand{\realdoublespacing}{%
  \let\CS=\@currsize\renewcommand{\baselinestretch}{2}\CS}
\newcommand{\listspacing}[1]{\singlespacing} % extra arg is extensibility hook

\newcommand{\prelim@contents}[1]{}
 \DeclareOption{tocprelim}{%
  \renewcommand{\prelim@contents}[1]{\addcontentsline{toc}{section}{#1}}
}
\ProcessOptions


% Packages
\usepackage[toc]{appendix}
\usepackage{fancyhdr}
\usepackage[letterpaper,left=1.5in,right=1.0in,top=1.0in,bottom=1.25in]{geometry}


% Page 
\PassOptionsToClass{\@fontsize}{report}%
\pagestyle{fancyplain}
\fancyhf{}
\rhead{ \fancyplain{\thepage}{\thepage} }

% Heading control
\renewcommand{\headrulewidth}{0pt} % remove lines as well
\renewcommand{\footrulewidth}{0pt}
\renewcommand{\prelim@contents}[1]{\addcontentsline{toc}{section}{#1}}
\renewcommand{\topfraction}{0.85}
\renewcommand{\textfraction}{0.1}
\renewcommand{\floatpagefraction}{0.75}
\renewcommand{\footnotesize}{\normalsize}
\def\@makechapterhead#1{%
  {\parindent \z@ \raggedright \normalfont
    \ifnum \c@secnumdepth >\m@ne
        \LARGE\bfseries \@chapapp\space \thechapter
        \par\nobreak
        \vskip 0\p@
    \fi
    \interlinepenalty\@M
    \Large \bfseries #1\par\nobreak
    \vskip 0\p@
  }}
\def\@makeschapterhead#1{%
  {
    \interlinepenalty\@M
    \LARGE \bfseries #1\par\nobreak
    \vskip 15\p@
  }}


% Chapter control
\renewcommand{\@makefntext}[1]{%
  \parindent 1em%
  \noindent
  \singlespacing
  \hbox to 1.8em{\hss\@makefnmark}#1
  \par\doublespacing
}

% Title page
\renewenvironment{titlepage}{%
  \typeout{Title page}
  \pagestyle{empty}
  \doublespacing
}{%
  \doublespacing\clearpage
}
\renewcommand{\maketitle}{%
  %\titlesize{\fontsize{12}{24}\selectfont}%
  \begin{titlepage}%
    \realdoublespacing
    \begin{center}%
      \hrule height 0.5in width 0pt % force vertical space at top
      University of Nevada, Reno
    \end{center}
    \vfill
    \begin{center}%
      \textbf{\@titlesize \expandafter\uppercase\expandafter{\@title}}
    \end{center}
    \vfill
    \begin{center}
      \@subtitlesize
      A \ifphd Dissertation \else Thesis \fi Submitted in Partial Fulfillment \\
      of the Requirements for the Degree of 
      \ifphd Doctor of Philosophy \else Master of \@masterstype \fi ~in \\
      \@degreefield
    \end{center}
    \vfill
    \begin{center}
      \@authorsize
      by \\
      \@author \\
      \@advisor ~/ \ifphd Dissertation \else Thesis \fi Advisor \\
      %\@conferralmonth\ \@conferralyear
    \end{center}
  \end{titlepage}
}


  
% Allow for the creation of the committee approval page
\newenvironment{committeepage}{%
  \typeout{Committee Approval page}
  \pagestyle{empty}
  \doublespacing
}{%
  \doublespacing\clearpage
}

% Larger spacing for the distance between committee member names
\newcommand{\doublespacingplus}{%
  \let\CS=\@currsize\renewcommand{\baselinestretch}{2.5}\CS}
  
% Set the text size for the header of the committee approval page
\newcommand{\@gradheadersize}{\fontsize{16}{24}\selectfont}
\newcommand{\gradheadersize}[1]{\renewcommand{\@gradheadersize}{#1}}

% Create the committee page
\newcommand{\makecommitteepage}{%
    \begin{committeepage}
        \noindent
        \makebox[\textwidth][c]{\textbf{\@gradheadersize THE GRADUATE SCHOOL}}
        \vfill
        \begin{center}%
            \singlespacing
            We recommend that the \ifphd dissertation \else thesis \fi prepared \\
            under our supervision by \\
            \doublespacing
        \end{center}
        \vfill
        \begin{center}%
            \textbf{\expandafter\uppercase\expandafter{\@author}} \\
        \end{center}
        \vfill
        \begin{center}%
            entitled
        \end{center}
        \vfill
        \begin{center}%
            \textbf{\expandafter\uppercase\expandafter{\@title}}
        \end{center}
        \vfill
        \begin{center}%
            \singlespacing
            be accepted in partial fulfillment of the \\
            requirements for the degree of \\
            \ifphd DOCTOR OF PHILOSOPHY \else MASTER OF 
            \expandafter\uppercase\expandafter{\@masterstype} \fi
            \doublespacing
        \end{center}
        \vfill
        \begin{center}%
            \doublespacingplus
            \listofmembers
            \doublespacing
        \end{center}
        \vfill
        \begin{center}%
            \@conferralmonth\ \@conferralyear
        \end{center}
    \end{committeepage}
}


\newenvironment{copyrightpage}{%
  \typeout{Copyleft page}
  \pagestyle{empty}
  \null\vfill
  \begin{center}
}{%
  \end{center}
  \vfill\null
  \newpage
  \clearpage%
}

\newcommand{\makecopyright}{%
  \begin{copyrightpage}
  \copyright\ \@copyrightyear{} \@copyrightholder \\
  ALL RIGHTS RESERVED
  \end{copyrightpage}%
}
\newcommand{\@copyrightholder}{\@author}
\newcommand{\@copyrightyear}{\@conferralyear}

\newcommand{\copyrightholder}[1]{\renewcommand{\@copyrightholder}{#1}}
\newcommand{\copyrightyear}[1]{\renewcommand{\@copyrightyear}{#1}}

\newcommand{\makepublicdomain}{%
  \begin{copyrightpage}
  This document is in the public domain.
  \end{copyrightpage}%
}


\newenvironment{acknowledgements}{%
  \prelim@contents{Acknowledgements}
  \chapter*{Acknowledgements}
}{%
  \cleardoublepage
}

\renewenvironment{abstract}{
  \doublespacing
  \pagenumbering{roman}
  \prelim@contents{Abstract}
  \chapter*{Abstract}
}{%
  \clearpage
}

% Environment for the main body of the thesis/dissertation
\newenvironment{manuscript}{
    \doublespacing \setcounter{page}{1} \pagenumbering{arabic}
    \addtolength{\parskip}{0.5\baselineskip}
}{
}

% Set the outpt style of the bibliography
\bibliographystyle{plain}

\renewenvironment{thebibliography}[1]
     {\chapter*{\bibname
        \@mkboth{\uppercase{\bibname}}{\uppercase{\bibname}}}%
      \list{\@biblabel{\arabic{enumiv}}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
    \itemsep=12pt
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\arabic{enumiv}}}%
      \sloppy\clubpenalty4000\widowpenalty4000%
      \sfcode`\.=\@m\relax
      \addcontentsline{toc}{chapter}{\bibname}
      \listspacing{bib}}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \doublespacing
      \endlist}

% Content formatting
\newcommand{\contentspage}{%
  \listspacing{toc}
  \label{toc}
  \prelim@contents{\contentsname}
  \tableofcontents
  \doublespacing
  \clearpage
}

\newcommand{\tablelistpage}{%
  \listspacing{tab}
  \label{lot}
  \prelim@contents{\listtablename}
  \listoftables
  \doublespacing
  \clearpage
}

\newcommand{\figurelistpage}{%
  \listspacing{fig}
  \label{lof}
  \prelim@contents{\listfigurename}
  \listoffigures
  \doublespacing
  \clearpage
}



\endinput
%%

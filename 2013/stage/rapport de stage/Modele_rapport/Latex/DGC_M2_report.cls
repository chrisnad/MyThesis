%% ====================================================================
%%  `LMTreport.cls' - a class for reports at LMT-Cachan
%%  -------------------------------------------------------------------
%%  David Néron
%%  neron@lmt.ens-cachan.fr
%% ====================================================================
\def\fileversion{1.3}
\def\filedate{2009/06/05}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{LMTreport}[\filedate\space v\fileversion\space-
                          a class for reports at LMT-Cachan]
\LoadClass[12pt]{book}

% ---------------------------------------> Chargement des packages <---

\RequirePackage{geometry}          % Régler la taille de la page
\RequirePackage{amsmath,amssymb,bm}% Étendre les fonctions maths
\RequirePackage{fancybox}          % Faire de jolies boîtes
\RequirePackage{calc}              % Faire des calculs
\RequirePackage{xspace}            % Gestion automatique des espaces
\RequirePackage{fancyhdr}          % Configurer headers/footers
\RequirePackage[french]{minitoc}   % Minitocs en début de chapitre
\RequirePackage{geometry}          % Géométrie de la page
\RequirePackage[centerlast,bf]{caption2} % Captions
\RequirePackage[center]{subfigure} % Inclusion de subfigures
\RequirePackage{array}             % Gestion des tableaux
\RequirePackage{tocbibind}         % ToC, ToF, etc. dans la ToC
\RequirePackage[frenchb]{babel}    % Règles typo françaises
\RequirePackage{ifpdf}             % Tester si génération d'un pdf

\usepackage{subfigure}

% ------------------------------------------> Géométrie de la page <---

\geometry{a4paper,
          twoside,
          body={15cm,22cm},
          lmargin=2.5cm,
          rmargin=2.5cm,
          footskip=1.5cm,
          top=3cm,
          headsep=1cm,
          head=15pt}

% ---------------------------------------------> LaTeX ou pdfLaTeX <---

\ifpdf
   \RequirePackage[pdftex]{graphicx,color}
   \RequirePackage[pdftex]{hyperref}
   \geometry{pdftex}
\else
   \RequirePackage[dvips]{graphicx,color}
   \RequirePackage[ps2pdf]{hyperref}
   \geometry{dvips}
\fi
	       
% --------------------------------------------> Réglage d'hyperref <---

\hypersetup{a4paper=true,
            pdfstartview=FitH,
	        colorlinks=false,
	        breaklinks=true,
	        pdfborder=0 0 0}

% ------------------------------------------------> toc et minitoc <---

% Réglage de la toc
\setcounter{tocdepth}{2}       % profondeur de la toc
\setcounter{secnumdepth}{3}    % profondeur de la numérotation 

% Réglage de la minitoc
\setcounter{minitocdepth}{2}
\setlength{\mtcindent}{15pt}                % indentation des minitocs
\renewcommand{\mtcfont}{\small\rm}          % fonte des minitocs
\renewcommand{\mtcSfont}{\small\bfseries}   % fonte des minitocs
\newcommand{\mtoc}{\minitoc\thispagestyle{plain}\newpage}

% Pour gérer des interférences entre minitoc et hyperref
\let\old@tableofcontents\tableofcontents
\def\tableofcontents{\old@tableofcontents\mtcaddchapter}
\let\old@listoffigures\listoffigures
\def\listoffigures{\old@listoffigures\mtcaddchapter}
\let\old@listoftables\listoftables
\def\listoftables{\old@listoftables\mtcaddchapter}

% --------------------------------------------> Quelques commandes <---

% Pour laisser un espace vertical entre deux paragraphes
\newcommand{\Par}{\bigskip\par}

% Pour que la page paire, ajoutée afin de commencer les chapitres
% sur une page impaire, soit sans header/footer
\newcommand{\clearemptydoublepage}{\newpage{\pagestyle{empty}%
                                   \cleardoublepage}}

% Des remerciements
\newenvironment{merci}
   {%
     \clearemptydoublepage
     \thispagestyle{empty}
     \null\vfill}
   {%
     \vfill\clearemptydoublepage}

% Des paragraphes encadrés
\newenvironment{parboxed}
   {%
     \setlength{\fboxsep}{10pt}
     \begin{Sbox}%
     \begin{minipage}{\linewidth-2\fboxrule-2\fboxsep}
     \setlength{\abovedisplayskip}{0pt}
     \setlength{\belowdisplayskip}{0pt}}
   {%
     \end{minipage}%
     \end{Sbox}
     \[\fbox{\TheSbox}\]\aftergroup\ignorespaces}

% Un résumé du chapitre
\newenvironment{resumechap}[1][]
   {%
     \vspace{2cm}
     \begin{center}
     \begin{minipage}{.7\linewidth}
     \centering\itshape}
   {%
     \end{minipage}
     \end{center}
     \vspace{2cm}}

% -----------------------------> Configuration des headers/footers <---

% Titre de bas de page
\def\lmt@foottitle{}
\newcommand{\foottitle}[1]{\def\lmt@foottitle{#1}}

% Pour une page courante : avec header et footer
\fancyhf{}
\fancyfoot[C]{\small\itshape\lmt@foottitle}
\fancyhead[RO,LE]{\small\bfseries\thepage}
\fancyhead[LO]{\nouppercase{\small\bfseries\rightmark}}
\fancyhead[RE]{\nouppercase{\small\bfseries\leftmark}}
\renewcommand{\footrulewidth}{.0pt}
\renewcommand{\headrulewidth}{.5pt}

% Pour une page de début de chapitre (ie de style plain) :
% avec seulement footer
\fancypagestyle{plain}{
\fancyhf{}
\fancyfoot[C]{\small\itshape\lmt@foottitle}
\renewcommand{\footrulewidth}{.0pt}
\renewcommand{\headrulewidth}{.0pt}}

\pagestyle{fancy}

\renewcommand{\chaptermark}[1]{\markboth{#1}{}}
\renewcommand{\sectionmark}[1]{\markright{#1}}

% -----------------------------------> Configuration des chapitres <---

% Numérotation moins chargée
\renewcommand{\thesection}{\@arabic\c@section}

% Style de chapitre, esthétique pompée dans thesetb.cls
\newfont{\chapfonte}{cmr10 scaled 9000}
\def\haut{20mm}
\def\decal{1.75ex}
\newlength{\lgnamechap}
\newlength{\lgnumchap}
\newlength{\lgall}    

% Redéfinition de \chapter
\def\@makechapterhead#1{%
  \vspace*{50\p@}%
  {\parindent \z@ \raggedleft \normalfont
    \ifnum \c@secnumdepth >\m@ne
      \if@mainmatter
        \huge\bfseries \@chapapp\space \thechapter
        \par\nobreak
        \vskip 20\p@
      \fi
    \fi
    \interlinepenalty\@M
    \Huge \bfseries #1\par\nobreak
    \vskip 40\p@
  }}

% Redéfinition de \chapter*
\def\@makeschapterhead#1{%
  \vspace*{50\p@}%
  {\parindent \z@ \raggedleft
    \normalfont
    \interlinepenalty\@M
    \Huge \bfseries  #1\par\nobreak
    \vskip 40\p@
  }\markboth{#1}{#1}}

% Chapitre sans numéro mais qui apparait dans la ToC
\newcommand{\chapterwithnonb}[1]{\toc@chapter{#1}\mtcaddchapter}

% ------------------------------> Gestion des flottants/références <---

% Noms des flottants
\newcommand{\Figure}{\textsc{Fig.}}
\newcommand{\Tableau}{\textsc{Tab.}}
\newcommand{\Algorithme}{\textsc{Alg.}}
\newcommand{\Chapitre}{\textsc{Chap.}}

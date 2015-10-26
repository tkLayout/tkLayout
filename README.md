


<!DOCTYPE html>
<html lang="en" class="">
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# object: http://ogp.me/ns/object# article: http://ogp.me/ns/article# profile: http://ogp.me/ns/profile#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta http-equiv="Content-Language" content="en">
    <meta name="viewport" content="width=1020">
    
    
    <title>tkLayout/README.md at devDelphes · drasal/tkLayout</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub">
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub">
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-114.png">
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114.png">
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-144.png">
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144.png">
    <meta property="fb:app_id" content="1401488693436528">

      <meta content="@github" name="twitter:site" /><meta content="summary" name="twitter:card" /><meta content="drasal/tkLayout" name="twitter:title" /><meta content="Contribute to tkLayout development by creating an account on GitHub." name="twitter:description" /><meta content="https://avatars0.githubusercontent.com/u/11534552?v=3&amp;s=400" name="twitter:image:src" />
      <meta content="GitHub" property="og:site_name" /><meta content="object" property="og:type" /><meta content="https://avatars0.githubusercontent.com/u/11534552?v=3&amp;s=400" property="og:image" /><meta content="drasal/tkLayout" property="og:title" /><meta content="https://github.com/drasal/tkLayout" property="og:url" /><meta content="Contribute to tkLayout development by creating an account on GitHub." property="og:description" />
      <meta name="browser-stats-url" content="https://api.github.com/_private/browser/stats">
    <meta name="browser-errors-url" content="https://api.github.com/_private/browser/errors">
    <link rel="assets" href="https://assets-cdn.github.com/">
    <link rel="web-socket" href="wss://live.github.com/_sockets/MTE1MzQ1NTI6ZTE3ZTE0MjYyMzA5NjcyZjU2MTAzZDUxNmU4YjlkZWQ6ZWQ0MDhhOTJhNThmY2RiNTY1ZDlhM2JhZTU4NjM5Nzg3YTBlNzk2OTVjZDE5MmM0NTc2NDNlZDY1MzFlYjIwMg==--c12c2597a9223161733cc5c294b0458f7d634aa0">
    <meta name="pjax-timeout" content="1000">
    <link rel="sudo-modal" href="/sessions/sudo_modal">

    <meta name="msapplication-TileImage" content="/windows-tile.png">
    <meta name="msapplication-TileColor" content="#ffffff">
    <meta name="selected-link" value="repo_source" data-pjax-transient>

    <meta name="google-site-verification" content="KT5gs8h0wvaagLKAVWq8bbeNwnZZK1r1XQysX3xurLU">
    <meta name="google-analytics" content="UA-3769691-2">

<meta content="collector.githubapp.com" name="octolytics-host" /><meta content="github" name="octolytics-app-id" /><meta content="808D25D8:3707:5FF3A81:562E7BC3" name="octolytics-dimension-request_id" /><meta content="11534552" name="octolytics-actor-id" /><meta content="drasal" name="octolytics-actor-login" /><meta content="ba059c06d3aafbc1e6e4ff7a34ef390446115f046d052f2eb987dfdde7c3d9de" name="octolytics-actor-hash" />

<meta content="Rails, view, blob#show" data-pjax-transient="true" name="analytics-event" />


  <meta class="js-ga-set" name="dimension1" content="Logged In">
    <meta class="js-ga-set" name="dimension4" content="Current repo nav">




    <meta name="is-dotcom" content="true">
        <meta name="hostname" content="github.com">
    <meta name="user-login" content="drasal">

      <link rel="mask-icon" href="https://assets-cdn.github.com/pinned-octocat.svg" color="#4078c0">
      <link rel="icon" type="image/x-icon" href="https://assets-cdn.github.com/favicon.ico">

    <meta content="a3df7e52d189ad3d7c114bc257dbceced158785e" name="form-nonce" />

    <link crossorigin="anonymous" href="https://assets-cdn.github.com/assets/github-fe9a78fd76e46beaf95a4330cc4974f836bfe46d4e95396090e339a58144486c.css" media="all" rel="stylesheet" />
    <link crossorigin="anonymous" href="https://assets-cdn.github.com/assets/github2-ed11b2e059c035eff96c8e0662efb7bbcffa8acc4978b0b94a8830b79ddf1736.css" media="all" rel="stylesheet" />
    
    
    


    <meta http-equiv="x-pjax-version" content="1788654b38a6bff8bcaf885c578af31d">

      
  <meta name="description" content="Contribute to tkLayout development by creating an account on GitHub.">
  <meta name="go-import" content="github.com/drasal/tkLayout git https://github.com/drasal/tkLayout.git">

  <meta content="11534552" name="octolytics-dimension-user_id" /><meta content="drasal" name="octolytics-dimension-user_login" /><meta content="34731135" name="octolytics-dimension-repository_id" /><meta content="drasal/tkLayout" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="true" name="octolytics-dimension-repository_is_fork" /><meta content="31327636" name="octolytics-dimension-repository_parent_id" /><meta content="tkLayout/tkLayout" name="octolytics-dimension-repository_parent_nwo" /><meta content="31327636" name="octolytics-dimension-repository_network_root_id" /><meta content="tkLayout/tkLayout" name="octolytics-dimension-repository_network_root_nwo" />
  <link href="https://github.com/drasal/tkLayout/commits/devDelphes.atom" rel="alternate" title="Recent Commits to tkLayout:devDelphes" type="application/atom+xml">

  </head>


  <body class="logged_in   env-production linux vis-public fork page-blob">
    <a href="#start-of-content" tabindex="1" class="accessibility-aid js-skip-to-content">Skip to content</a>

    
    
    



      <div class="header header-logged-in true" role="banner">
  <div class="container clearfix">

    <a class="header-logo-invertocat" href="https://github.com/" data-hotkey="g d" aria-label="Homepage" data-ga-click="Header, go to dashboard, icon:logo">
  <span class="mega-octicon octicon-mark-github"></span>
</a>


      <div class="site-search repo-scope js-site-search" role="search">
          <!-- </textarea> --><!-- '"` --><form accept-charset="UTF-8" action="/drasal/tkLayout/search" class="js-site-search-form" data-global-search-url="/search" data-repo-search-url="/drasal/tkLayout/search" method="get"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /></div>
  <label class="js-chromeless-input-container form-control">
    <div class="scope-badge">This repository</div>
    <input type="text"
      class="js-site-search-focus js-site-search-field is-clearable chromeless-input"
      data-hotkey="s"
      name="q"
      placeholder="Search"
      aria-label="Search this repository"
      data-global-scope-placeholder="Search GitHub"
      data-repo-scope-placeholder="Search"
      tabindex="1"
      autocapitalize="off">
  </label>
</form>
      </div>

      <ul class="header-nav left" role="navigation">
        <li class="header-nav-item">
          <a href="/pulls" class="js-selected-navigation-item header-nav-link" data-ga-click="Header, click, Nav menu - item:pulls context:user" data-hotkey="g p" data-selected-links="/pulls /pulls/assigned /pulls/mentioned /pulls">
            Pull requests
</a>        </li>
        <li class="header-nav-item">
          <a href="/issues" class="js-selected-navigation-item header-nav-link" data-ga-click="Header, click, Nav menu - item:issues context:user" data-hotkey="g i" data-selected-links="/issues /issues/assigned /issues/mentioned /issues">
            Issues
</a>        </li>
          <li class="header-nav-item">
            <a class="header-nav-link" href="https://gist.github.com/" data-ga-click="Header, go to gist, text:gist">Gist</a>
          </li>
      </ul>

    
<ul class="header-nav user-nav right" id="user-links">
  <li class="header-nav-item">
      <span class="js-socket-channel js-updatable-content"
        data-channel="notification-changed:drasal"
        data-url="/notifications/header">
      <a href="/notifications" aria-label="You have unread notifications" class="header-nav-link notification-indicator tooltipped tooltipped-s" data-ga-click="Header, go to notifications, icon:unread" data-hotkey="g n">
          <span class="mail-status unread"></span>
          <span class="octicon octicon-bell"></span>
</a>  </span>

  </li>

  <li class="header-nav-item dropdown js-menu-container">
    <a class="header-nav-link tooltipped tooltipped-s js-menu-target" href="/new"
       aria-label="Create new…"
       data-ga-click="Header, create new, icon:add">
      <span class="octicon octicon-plus left"></span>
      <span class="dropdown-caret"></span>
    </a>

    <div class="dropdown-menu-content js-menu-content">
      <ul class="dropdown-menu dropdown-menu-sw">
        
<a class="dropdown-item" href="/new" data-ga-click="Header, create new repository">
  New repository
</a>


  <a class="dropdown-item" href="/organizations/new" data-ga-click="Header, create new organization">
    New organization
  </a>



  <div class="dropdown-divider"></div>
  <div class="dropdown-header">
    <span title="drasal/tkLayout">This repository</span>
  </div>
    <a class="dropdown-item" href="/drasal/tkLayout/settings/collaboration" data-ga-click="Header, create new collaborator">
      New collaborator
    </a>

      </ul>
    </div>
  </li>

  <li class="header-nav-item dropdown js-menu-container">
    <a class="header-nav-link name tooltipped tooltipped-s js-menu-target" href="/drasal"
       aria-label="View profile and more"
       data-ga-click="Header, show menu, icon:avatar">
      <img alt="@drasal" class="avatar" height="20" src="https://avatars1.githubusercontent.com/u/11534552?v=3&amp;s=40" width="20" />
      <span class="dropdown-caret"></span>
    </a>

    <div class="dropdown-menu-content js-menu-content">
      <div class="dropdown-menu  dropdown-menu-sw">
        <div class=" dropdown-header header-nav-current-user css-truncate">
            Signed in as <strong class="css-truncate-target">drasal</strong>

        </div>


        <div class="dropdown-divider"></div>

          <a class="dropdown-item" href="/drasal" data-ga-click="Header, go to profile, text:your profile">
            Your profile
          </a>
        <a class="dropdown-item" href="/stars" data-ga-click="Header, go to starred repos, text:your stars">
          Your stars
        </a>
        <a class="dropdown-item" href="/explore" data-ga-click="Header, go to explore, text:explore">
          Explore
        </a>
          <a class="dropdown-item" href="/integrations" data-ga-click="Header, go to integrations, text:integrations">
            Integrations
          </a>
        <a class="dropdown-item" href="https://help.github.com" data-ga-click="Header, go to help, text:help">
          Help
        </a>

          <div class="dropdown-divider"></div>

          <a class="dropdown-item" href="/settings/profile" data-ga-click="Header, go to settings, icon:settings">
            Settings
          </a>

          <!-- </textarea> --><!-- '"` --><form accept-charset="UTF-8" action="/logout" class="logout-form" data-form-nonce="a3df7e52d189ad3d7c114bc257dbceced158785e" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="iuqhjQZcHwNCSI3UAPKe6vbtJ9piu3qS9+f/g5ek5DLhAS71h+txlExrsu9pAP317Etb6iYljwBPdU3wtkPq7w==" /></div>
            <button class="dropdown-item dropdown-signout" data-ga-click="Header, sign out, icon:logout">
              Sign out
            </button>
</form>
      </div>
    </div>
  </li>
</ul>


    
  </div>
</div>

      

      


    <div id="start-of-content" class="accessibility-aid"></div>

    <div id="js-flash-container">
</div>


    <div role="main" class="main-content">
        <div itemscope itemtype="http://schema.org/WebPage">
    <div class="pagehead repohead instapaper_ignore readability-menu">

      <div class="container">

        <div class="clearfix">
          

<ul class="pagehead-actions">

  <li>
      <!-- </textarea> --><!-- '"` --><form accept-charset="UTF-8" action="/notifications/subscribe" class="js-social-container" data-autosubmit="true" data-form-nonce="a3df7e52d189ad3d7c114bc257dbceced158785e" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="/zdWJOG0cmlSPzcgUD4wv8xVAPhCDa9gesFD+PighI2HHdxvRZnL2iksHEmwiAWyXsSDnkCstKsxmipJWuIHhw==" /></div>    <input id="repository_id" name="repository_id" type="hidden" value="34731135" />

      <div class="select-menu js-menu-container js-select-menu">
        <a href="/drasal/tkLayout/subscription"
          class="btn btn-sm btn-with-count select-menu-button js-menu-target" role="button" tabindex="0" aria-haspopup="true"
          data-ga-click="Repository, click Watch settings, action:blob#show">
          <span class="js-select-button">
            <span class="octicon octicon-eye"></span>
            Unwatch
          </span>
        </a>
        <a class="social-count js-social-count" href="/drasal/tkLayout/watchers">
          1
        </a>

        <div class="select-menu-modal-holder">
          <div class="select-menu-modal subscription-menu-modal js-menu-content" aria-hidden="true">
            <div class="select-menu-header">
              <span class="select-menu-title">Notifications</span>
              <span class="octicon octicon-x js-menu-close" role="button" aria-label="Close"></span>
            </div>

            <div class="select-menu-list js-navigation-container" role="menu">

              <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
                <span class="select-menu-item-icon octicon octicon-check"></span>
                <div class="select-menu-item-text">
                  <input id="do_included" name="do" type="radio" value="included" />
                  <span class="select-menu-item-heading">Not watching</span>
                  <span class="description">Be notified when participating or @mentioned.</span>
                  <span class="js-select-button-text hidden-select-button-text">
                    <span class="octicon octicon-eye"></span>
                    Watch
                  </span>
                </div>
              </div>

              <div class="select-menu-item js-navigation-item selected" role="menuitem" tabindex="0">
                <span class="select-menu-item-icon octicon octicon octicon-check"></span>
                <div class="select-menu-item-text">
                  <input checked="checked" id="do_subscribed" name="do" type="radio" value="subscribed" />
                  <span class="select-menu-item-heading">Watching</span>
                  <span class="description">Be notified of all conversations.</span>
                  <span class="js-select-button-text hidden-select-button-text">
                    <span class="octicon octicon-eye"></span>
                    Unwatch
                  </span>
                </div>
              </div>

              <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
                <span class="select-menu-item-icon octicon octicon-check"></span>
                <div class="select-menu-item-text">
                  <input id="do_ignore" name="do" type="radio" value="ignore" />
                  <span class="select-menu-item-heading">Ignoring</span>
                  <span class="description">Never be notified.</span>
                  <span class="js-select-button-text hidden-select-button-text">
                    <span class="octicon octicon-mute"></span>
                    Stop ignoring
                  </span>
                </div>
              </div>

            </div>

          </div>
        </div>
      </div>
</form>
  </li>

  <li>
    
  <div class="js-toggler-container js-social-container starring-container ">

    <!-- </textarea> --><!-- '"` --><form accept-charset="UTF-8" action="/drasal/tkLayout/unstar" class="js-toggler-form starred js-unstar-button" data-form-nonce="a3df7e52d189ad3d7c114bc257dbceced158785e" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="gPlcMGHR0ZJJY4bbSLJXpB3U15XM1zsTBbHj7nNE8QMXNpmaYyl2a/QretBUn6MRJAqqtn2m3bsT99PevCWp4w==" /></div>
      <button
        class="btn btn-sm btn-with-count js-toggler-target"
        aria-label="Unstar this repository" title="Unstar drasal/tkLayout"
        data-ga-click="Repository, click unstar button, action:blob#show; text:Unstar">
        <span class="octicon octicon-star"></span>
        Unstar
      </button>
        <a class="social-count js-social-count" href="/drasal/tkLayout/stargazers">
          0
        </a>
</form>
    <!-- </textarea> --><!-- '"` --><form accept-charset="UTF-8" action="/drasal/tkLayout/star" class="js-toggler-form unstarred js-star-button" data-form-nonce="a3df7e52d189ad3d7c114bc257dbceced158785e" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="iDioXEyhktnwjSbH0iKa76c44ghoKTUHbrtbd1TbgQvkN6Ptliai+UYWJmjVZfQOwwQDxlLBunpxbhv8M41jeQ==" /></div>
      <button
        class="btn btn-sm btn-with-count js-toggler-target"
        aria-label="Star this repository" title="Star drasal/tkLayout"
        data-ga-click="Repository, click star button, action:blob#show; text:Star">
        <span class="octicon octicon-star"></span>
        Star
      </button>
        <a class="social-count js-social-count" href="/drasal/tkLayout/stargazers">
          0
        </a>
</form>  </div>

  </li>

  <li>
          <a href="#fork-destination-box" class="btn btn-sm btn-with-count"
              title="Fork your own copy of drasal/tkLayout to your account"
              aria-label="Fork your own copy of drasal/tkLayout to your account"
              rel="facebox"
              data-ga-click="Repository, show fork modal, action:blob#show; text:Fork">
            <span class="octicon octicon-repo-forked"></span>
            Fork
          </a>

          <div id="fork-destination-box" style="display: none;">
            <h2 class="facebox-header" data-facebox-id="facebox-header">Where should we fork this repository?</h2>
            <include-fragment src=""
                class="js-fork-select-fragment fork-select-fragment"
                data-url="/drasal/tkLayout/fork?fragment=1">
              <img alt="Loading" height="64" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-128.gif" width="64" />
            </include-fragment>
          </div>

    <a href="/drasal/tkLayout/network" class="social-count">
      10
    </a>
  </li>
</ul>

          <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public ">
  <span class="mega-octicon octicon-repo-forked"></span>
  <span class="author"><a href="/drasal" class="url fn" itemprop="url" rel="author"><span itemprop="title">drasal</span></a></span><!--
--><span class="path-divider">/</span><!--
--><strong><a href="/drasal/tkLayout" data-pjax="#js-repo-pjax-container">tkLayout</a></strong>

  <span class="page-context-loader">
    <img alt="" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
  </span>

    <span class="fork-flag">
      <span class="text">forked from <a href="/tkLayout/tkLayout">tkLayout/tkLayout</a></span>
    </span>
</h1>

        </div>
      </div>
    </div>

    <div class="container">
      <div class="repository-with-sidebar repo-container new-discussion-timeline ">
        <div class="repository-sidebar clearfix">
          
<nav class="sunken-menu repo-nav js-repo-nav js-sidenav-container-pjax js-octicon-loaders"
     role="navigation"
     data-pjax="#js-repo-pjax-container"
     data-issue-count-url="/drasal/tkLayout/issues/counts">
  <ul class="sunken-menu-group">
    <li class="tooltipped tooltipped-w" aria-label="Code">
      <a href="/drasal/tkLayout/tree/devDelphes" aria-label="Code" aria-selected="true" class="js-selected-navigation-item selected sunken-menu-item" data-hotkey="g c" data-selected-links="repo_source repo_downloads repo_commits repo_releases repo_tags repo_branches /drasal/tkLayout/tree/devDelphes">
        <span class="octicon octicon-code"></span> <span class="full-word">Code</span>
        <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>    </li>


    <li class="tooltipped tooltipped-w" aria-label="Pull requests">
      <a href="/drasal/tkLayout/pulls" aria-label="Pull requests" class="js-selected-navigation-item sunken-menu-item" data-hotkey="g p" data-selected-links="repo_pulls /drasal/tkLayout/pulls">
          <span class="octicon octicon-git-pull-request"></span> <span class="full-word">Pull requests</span>
          <span class="js-pull-replace-counter"></span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>    </li>

      <li class="tooltipped tooltipped-w" aria-label="Wiki">
        <a href="/drasal/tkLayout/wiki" aria-label="Wiki" class="js-selected-navigation-item sunken-menu-item" data-hotkey="g w" data-selected-links="repo_wiki /drasal/tkLayout/wiki">
          <span class="octicon octicon-book"></span> <span class="full-word">Wiki</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>
  </ul>
  <div class="sunken-menu-separator"></div>
  <ul class="sunken-menu-group">

    <li class="tooltipped tooltipped-w" aria-label="Pulse">
      <a href="/drasal/tkLayout/pulse" aria-label="Pulse" class="js-selected-navigation-item sunken-menu-item" data-selected-links="pulse /drasal/tkLayout/pulse">
        <span class="octicon octicon-pulse"></span> <span class="full-word">Pulse</span>
        <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>    </li>

    <li class="tooltipped tooltipped-w" aria-label="Graphs">
      <a href="/drasal/tkLayout/graphs" aria-label="Graphs" class="js-selected-navigation-item sunken-menu-item" data-selected-links="repo_graphs repo_contributors /drasal/tkLayout/graphs">
        <span class="octicon octicon-graph"></span> <span class="full-word">Graphs</span>
        <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>    </li>
  </ul>


    <div class="sunken-menu-separator"></div>
    <ul class="sunken-menu-group">
      <li class="tooltipped tooltipped-w" aria-label="Settings">
        <a href="/drasal/tkLayout/settings" aria-label="Settings" class="js-selected-navigation-item sunken-menu-item" data-selected-links="repo_settings repo_branch_settings hooks /drasal/tkLayout/settings">
          <span class="octicon octicon-gear"></span> <span class="full-word">Settings</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>
    </ul>
</nav>

            <div class="only-with-full-nav">
                
<div class="js-clone-url clone-url open"
  data-protocol-type="http">
  <h3 class="text-small text-muted"><span class="text-emphasized">HTTPS</span> clone URL</h3>
  <div class="input-group js-zeroclipboard-container">
    <input type="text" class="input-mini text-small text-muted input-monospace js-url-field js-zeroclipboard-target"
           value="https://github.com/drasal/tkLayout.git" readonly="readonly" aria-label="HTTPS clone URL">
    <span class="input-group-button">
      <button aria-label="Copy to clipboard" class="js-zeroclipboard btn btn-sm zeroclipboard-button tooltipped tooltipped-s" data-copied-hint="Copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>

  
<div class="js-clone-url clone-url "
  data-protocol-type="ssh">
  <h3 class="text-small text-muted"><span class="text-emphasized">SSH</span> clone URL</h3>
  <div class="input-group js-zeroclipboard-container">
    <input type="text" class="input-mini text-small text-muted input-monospace js-url-field js-zeroclipboard-target"
           value="git@github.com:drasal/tkLayout.git" readonly="readonly" aria-label="SSH clone URL">
    <span class="input-group-button">
      <button aria-label="Copy to clipboard" class="js-zeroclipboard btn btn-sm zeroclipboard-button tooltipped tooltipped-s" data-copied-hint="Copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>

  
<div class="js-clone-url clone-url "
  data-protocol-type="subversion">
  <h3 class="text-small text-muted"><span class="text-emphasized">Subversion</span> checkout URL</h3>
  <div class="input-group js-zeroclipboard-container">
    <input type="text" class="input-mini text-small text-muted input-monospace js-url-field js-zeroclipboard-target"
           value="https://github.com/drasal/tkLayout" readonly="readonly" aria-label="Subversion checkout URL">
    <span class="input-group-button">
      <button aria-label="Copy to clipboard" class="js-zeroclipboard btn btn-sm zeroclipboard-button tooltipped tooltipped-s" data-copied-hint="Copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>



<div class="clone-options text-small text-muted">You can clone with
  <!-- </textarea> --><!-- '"` --><form accept-charset="UTF-8" action="/users/set_protocol?protocol_selector=http&amp;protocol_type=push" class="inline-form js-clone-selector-form is-enabled" data-form-nonce="a3df7e52d189ad3d7c114bc257dbceced158785e" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="RD8BK1xZ+CCHKaAlK0hBYYLlDtZyILqziZQDiGPMLWBi3NSmGH3ZPK8srKHtsRCNPQLWmluInH5QCIALlktSEw==" /></div><button class="btn-link js-clone-selector" data-protocol="http" type="submit">HTTPS</button></form>, <!-- </textarea> --><!-- '"` --><form accept-charset="UTF-8" action="/users/set_protocol?protocol_selector=ssh&amp;protocol_type=push" class="inline-form js-clone-selector-form is-enabled" data-form-nonce="a3df7e52d189ad3d7c114bc257dbceced158785e" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="J3W2j7zsXdq8ikSjImOo8lfHU6BOAe8thOQwF50AOk1sXl0k+YHHDjQ2O/ehfkrQ7p7aFInZJu5bQWm+kh3deQ==" /></div><button class="btn-link js-clone-selector" data-protocol="ssh" type="submit">SSH</button></form>, or <!-- </textarea> --><!-- '"` --><form accept-charset="UTF-8" action="/users/set_protocol?protocol_selector=subversion&amp;protocol_type=push" class="inline-form js-clone-selector-form is-enabled" data-form-nonce="a3df7e52d189ad3d7c114bc257dbceced158785e" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="Gou8pCES1iaI0EHNLTYDS7TlMrMA8w5R68xKqYaQGduXz+xIQz+kH9cfLyfaGH8VUjLbgBNfc6iek0e38eCdAA==" /></div><button class="btn-link js-clone-selector" data-protocol="subversion" type="submit">Subversion</button></form>.
  <a href="https://help.github.com/articles/which-remote-url-should-i-use" class="help tooltipped tooltipped-n" aria-label="Get help on which URL is right for you.">
    <span class="octicon octicon-question"></span>
  </a>
</div>

              <a href="/drasal/tkLayout/archive/devDelphes.zip"
                 class="btn btn-sm sidebar-button"
                 aria-label="Download the contents of drasal/tkLayout as a zip file"
                 title="Download the contents of drasal/tkLayout as a zip file"
                 rel="nofollow">
                <span class="octicon octicon-cloud-download"></span>
                Download ZIP
              </a>
            </div>
        </div>
        <div id="js-repo-pjax-container" class="repository-content context-loader-container" data-pjax-container>

          

<a href="/drasal/tkLayout/blob/3551c5cd56c25c263c5a7a4176ab952eae69363e/README.md" class="hidden js-permalink-shortcut" data-hotkey="y">Permalink</a>

<!-- blob contrib key: blob_contributors:v21:a3a4d5ade552c668b0946373e024b7f0 -->

  <div class="file-navigation js-zeroclipboard-container">
    
<div class="select-menu js-menu-container js-select-menu left">
  <button class="btn btn-sm select-menu-button js-menu-target css-truncate" data-hotkey="w"
    title="devDelphes"
    type="button" aria-label="Switch branches or tags" tabindex="0" aria-haspopup="true">
    <i>Branch:</i>
    <span class="js-select-button css-truncate-target">devDelphes</span>
  </button>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax aria-hidden="true">

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <span class="select-menu-title">Switch branches/tags</span>
        <span class="octicon octicon-x js-menu-close" role="button" aria-label="Close"></span>
      </div>

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" aria-label="Find or create a branch…" id="context-commitish-filter-field" class="js-filterable-field js-navigation-enable" placeholder="Find or create a branch…">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" data-filter-placeholder="Find or create a branch…" class="js-select-menu-tab" role="tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" data-filter-placeholder="Find a tag…" class="js-select-menu-tab" role="tab">Tags</a>
            </li>
          </ul>
        </div>
      </div>

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches" role="menu">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/CUPS2014/README.md"
               data-name="CUPS2014"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="CUPS2014">
                CUPS2014
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/cmssw22x/README.md"
               data-name="cmssw22x"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="cmssw22x">
                cmssw22x
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/coderev/README.md"
               data-name="coderev"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="coderev">
                coderev
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/cups2014.2/README.md"
               data-name="cups2014.2"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="cups2014.2">
                cups2014.2
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/deMaioMaterial/README.md"
               data-name="deMaioMaterial"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="deMaioMaterial">
                deMaioMaterial
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/dev/README.md"
               data-name="dev"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="dev">
                dev
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/devCMake/README.md"
               data-name="devCMake"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="devCMake">
                devCMake
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open selected"
               href="/drasal/tkLayout/blob/devDelphes/README.md"
               data-name="devDelphes"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="devDelphes">
                devDelphes
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/devFwd/README.md"
               data-name="devFwd"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="devFwd">
                devFwd
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/devPMom/README.md"
               data-name="devPMom"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="devPMom">
                devPMom
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/fcc/README.md"
               data-name="fcc"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="fcc">
                fcc
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/gbianchi/README.md"
               data-name="gbianchi"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="gbianchi">
                gbianchi
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/gbianchi2/README.md"
               data-name="gbianchi2"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="gbianchi2">
                gbianchi2
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/jelena/README.md"
               data-name="jelena"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="jelena">
                jelena
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/master/README.md"
               data-name="master"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="master">
                master
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/materialLayer/README.md"
               data-name="materialLayer"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="materialLayer">
                materialLayer
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/ryo/README.md"
               data-name="ryo"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="ryo">
                ryo
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/stable/README.md"
               data-name="stable"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="stable">
                stable
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/tiltmodgeom/README.md"
               data-name="tiltmodgeom"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="tiltmodgeom">
                tiltmodgeom
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/usherRefactoring/README.md"
               data-name="usherRefactoring"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="usherRefactoring">
                usherRefactoring
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/drasal/tkLayout/blob/zfit/README.md"
               data-name="zfit"
               data-skip-pjax="true"
               rel="nofollow">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <span class="select-menu-item-text css-truncate-target" title="zfit">
                zfit
              </span>
            </a>
        </div>

          <!-- </textarea> --><!-- '"` --><form accept-charset="UTF-8" action="/drasal/tkLayout/branches" class="js-create-branch select-menu-item select-menu-new-item-form js-navigation-item js-new-item-form" data-form-nonce="a3df7e52d189ad3d7c114bc257dbceced158785e" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="DWTc5tevP4et5AMYerNtaoFFFZVebe+8jVV1sbamKAUsNTcv0SQtsd+01BMWaBC/1BH1kTTyhX0qWqMGQZyhWg==" /></div>
            <span class="octicon octicon-git-branch select-menu-item-icon"></span>
            <div class="select-menu-item-text">
              <span class="select-menu-item-heading">Create branch: <span class="js-new-item-name"></span></span>
              <span class="description">from ‘devDelphes’</span>
            </div>
            <input type="hidden" name="name" id="name" class="js-new-item-value">
            <input type="hidden" name="branch" id="branch" value="devDelphes">
            <input type="hidden" name="path" id="path" value="README.md">
</form>
      </div>

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/drasal/tkLayout/tree/4.0/README.md"
                 data-name="4.0"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text css-truncate-target"
                 title="4.0">4.0</a>
            </div>
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/drasal/tkLayout/tree/3.614/README.md"
                 data-name="3.614"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text css-truncate-target"
                 title="3.614">3.614</a>
            </div>
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/drasal/tkLayout/tree/3.9_TP/README.md"
                 data-name="3.9_TP"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text css-truncate-target"
                 title="3.9_TP">3.9_TP</a>
            </div>
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/drasal/tkLayout/tree/3.8/README.md"
                 data-name="3.8"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text css-truncate-target"
                 title="3.8">3.8</a>
            </div>
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/drasal/tkLayout/tree/3.0/README.md"
                 data-name="3.0"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text css-truncate-target"
                 title="3.0">3.0</a>
            </div>
        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div>

    </div>
  </div>
</div>

    <div class="btn-group right">
      <a href="/drasal/tkLayout/find/devDelphes"
            class="js-show-file-finder btn btn-sm empty-icon tooltipped tooltipped-nw"
            data-pjax
            data-hotkey="t"
            aria-label="Quickly jump between files">
        <span class="octicon octicon-list-unordered"></span>
      </a>
      <button aria-label="Copy file path to clipboard" class="js-zeroclipboard btn btn-sm zeroclipboard-button tooltipped tooltipped-s" data-copied-hint="Copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </div>

    <div class="breadcrumb js-zeroclipboard-target">
      <span class="repo-root js-repo-root"><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/drasal/tkLayout/tree/devDelphes" class="" data-branch="devDelphes" data-pjax="true" itemscope="url"><span itemprop="title">tkLayout</span></a></span></span><span class="separator">/</span><strong class="final-path">README.md</strong>
    </div>
  </div>


  <div class="commit-tease">
      <span class="right">
        <a class="commit-tease-sha" href="/drasal/tkLayout/commit/ab9c053b7c97e06ca4dab3f0ce6a515b08cf8bce" data-pjax>
          ab9c053
        </a>
        <time datetime="2015-06-10T08:18:03Z" is="relative-time">Jun 10, 2015</time>
      </span>
      <div>
        <img alt="@drasal" class="avatar" height="20" src="https://avatars1.githubusercontent.com/u/11534552?v=3&amp;s=40" width="20" />
        <a href="/drasal" class="user-mention" rel="author">drasal</a>
          <a href="/drasal/tkLayout/commit/ab9c053b7c97e06ca4dab3f0ce6a515b08cf8bce" class="message" data-pjax="true" title="Added Cross-Platform compiling system - CMake, follow README.md (rearanged, added comments) ...">Added Cross-Platform compiling system - CMake, follow README.md (rear…</a>
      </div>

    <div class="commit-tease-contributors">
      <a class="muted-link contributors-toggle" href="#blob_contributors_box" rel="facebox">
        <strong>1</strong>
         contributor
      </a>
      
    </div>

    <div id="blob_contributors_box" style="display:none">
      <h2 class="facebox-header" data-facebox-id="facebox-header">Users who have contributed to this file</h2>
      <ul class="facebox-user-list" data-facebox-id="facebox-description">
          <li class="facebox-user-list-item">
            <img alt="@drasal" height="24" src="https://avatars3.githubusercontent.com/u/11534552?v=3&amp;s=48" width="24" />
            <a href="/drasal">drasal</a>
          </li>
      </ul>
    </div>
  </div>

<div class="file">
  <div class="file-header">
  <div class="file-actions">

    <div class="btn-group">
      <a href="/drasal/tkLayout/raw/devDelphes/README.md" class="btn btn-sm " id="raw-url">Raw</a>
        <a href="/drasal/tkLayout/blame/devDelphes/README.md" class="btn btn-sm js-update-url-with-hash">Blame</a>
      <a href="/drasal/tkLayout/commits/devDelphes/README.md" class="btn btn-sm " rel="nofollow">History</a>
    </div>


        <!-- </textarea> --><!-- '"` --><form accept-charset="UTF-8" action="/drasal/tkLayout/edit/devDelphes/README.md" class="inline-form js-update-url-with-hash" data-form-nonce="a3df7e52d189ad3d7c114bc257dbceced158785e" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="2P20ecHgVMqVXmOuCPjsut7XMvn1ddGXCK89WcMeg+ItrXdeiTS01UAMgJd9KT4u7AVTKBXpLvxmAQtWLk9OgQ==" /></div>
          <button class="octicon-btn tooltipped tooltipped-nw" type="submit"
            aria-label="Edit this file" data-hotkey="e" data-disable-with>
            <span class="octicon octicon-pencil"></span>
          </button>
</form>        <!-- </textarea> --><!-- '"` --><form accept-charset="UTF-8" action="/drasal/tkLayout/delete/devDelphes/README.md" class="inline-form" data-form-nonce="a3df7e52d189ad3d7c114bc257dbceced158785e" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="T3Mxom99r5lQ+iAMTOSlwGaL3JpECZFVxX2nixWPHzouUfaSU03gVAJzLFgndsbF2hR1VE1cOBt9N29AV0RkcQ==" /></div>
          <button class="octicon-btn octicon-btn-danger tooltipped tooltipped-nw" type="submit"
            aria-label="Delete this file" data-disable-with>
            <span class="octicon octicon-trashcan"></span>
          </button>
</form>  </div>

  <div class="file-info">
      69 lines (46 sloc)
      <span class="file-info-divider"></span>
    2.27 KB
  </div>
</div>

  
  <div id="readme" class="blob instapaper_body">
    <article class="markdown-body entry-content" itemprop="mainContentOfPage"><h1><a id="user-content-getting-the-code" class="anchor" href="#getting-the-code" aria-hidden="true"><span class="octicon octicon-link"></span></a>Getting the code</h1>

<p>git clone <a href="https://github.com/tkLayout/tkLayout.git">https://github.com/tkLayout/tkLayout.git</a>
   cd tkLayout</p>

<h1><a id="user-content-before-the-compilationrun" class="anchor" href="#before-the-compilationrun" aria-hidden="true"><span class="octicon octicon-link"></span></a>Before the compilation/run</h1>

<p>You need a working version of root active (root and root-config should be in your path) and a few
libraries are required:</p>

<ul>
<li>boost_filesystem</li>
<li>boost_regex</li>
<li>the root library set</li>
</ul>

<p>If you are on lxplus you just have to run a bash shell and then</p>

<pre><code> source setup_slc6.sh
</code></pre>

<h1><a id="user-content-compilation" class="anchor" href="#compilation" aria-hidden="true"><span class="octicon octicon-link"></span></a>Compilation</h1>

<p>Compilation using make:</p>

<pre><code>make
</code></pre>

<p>This will build the needed programs and put them in the ./bin directory.</p>

<p>Compilation using cmake (NEW):</p>

<pre><code>mkdir build     (all object files, help files, libs, ... will be kept here)
cd build
cmake ..        (generate makefile)
make install    (or write make all + make install)

make uninstall  (if cleaning needed)
rm *            (clean all content in build directory &amp; restart if needed)
</code></pre>

<p>This will build the needed programs, copy the executables into tkLayout/bin directory and create symlink in ${HOME}/bin directory.</p>

<h1><a id="user-content-install" class="anchor" href="#install" aria-hidden="true"><span class="octicon octicon-link"></span></a>Install</h1>

<p>If the make command runs properly you can install the program with the script</p>

<pre><code> make install
</code></pre>

<h2><a id="user-content-first-time-install" class="anchor" href="#first-time-install" aria-hidden="true"><span class="octicon octicon-link"></span></a>First-time install</h2>

<p>If this is the first time that you install tkLayout, a few questions will be asked and a configuration
fle will be created in $HOME/.tkgeometry:</p>

<ol>
<li>You will need to provide the destination directory for the www output of the program (proabably
something like /afs/cern.ch/user/y/yourname/www/layouts) this directory needs to be writable by you
and readable by the web server.
   The style files will be copied during the installation to that directory (for example
/afs/cern.ch/user/y/yourname/www/layouts/style ). If the style files get changed during the development
a new ./install.sh will be needed to get this propagated to the output.
   Alternatively, you can choose to make a symbolic link between the source directory style directory and
the layout directory. This avoids the need of repeating the ./install.sh at every update of the style files
but in order to do this the source directory should also be within the reach of the web server.</li>
<li>The set of momentum values to be used for the simulation</li>
<li>[...] t.b.d.</li>
</ol>

<h1><a id="user-content-update" class="anchor" href="#update" aria-hidden="true"><span class="octicon octicon-link"></span></a>Update</h1>

<p>To get the latest development version (usually UNSTABLE) you simply type</p>

<pre><code>git fetch
git chechout master
make
make install
</code></pre>
</article>
  </div>

</div>

<a href="#jump-to-line" rel="facebox[.linejump]" data-hotkey="l" style="display:none">Jump to Line</a>
<div id="jump-to-line" style="display:none">
  <!-- </textarea> --><!-- '"` --><form accept-charset="UTF-8" action="" class="js-jump-to-line-form" method="get"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /></div>
    <input class="linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" aria-label="Jump to line" autofocus>
    <button type="submit" class="btn">Go</button>
</form></div>

        </div>
      </div>
      <div class="modal-backdrop"></div>
    </div>
  </div>


    </div>

      <div class="container">
  <div class="site-footer" role="contentinfo">
    <ul class="site-footer-links right">
        <li><a href="https://status.github.com/" data-ga-click="Footer, go to status, text:status">Status</a></li>
      <li><a href="https://developer.github.com" data-ga-click="Footer, go to api, text:api">API</a></li>
      <li><a href="https://training.github.com" data-ga-click="Footer, go to training, text:training">Training</a></li>
      <li><a href="https://shop.github.com" data-ga-click="Footer, go to shop, text:shop">Shop</a></li>
        <li><a href="https://github.com/blog" data-ga-click="Footer, go to blog, text:blog">Blog</a></li>
        <li><a href="https://github.com/about" data-ga-click="Footer, go to about, text:about">About</a></li>
        <li><a href="https://github.com/pricing" data-ga-click="Footer, go to pricing, text:pricing">Pricing</a></li>

    </ul>

    <a href="https://github.com" aria-label="Homepage">
      <span class="mega-octicon octicon-mark-github" title="GitHub"></span>
</a>
    <ul class="site-footer-links">
      <li>&copy; 2015 <span title="0.12603s from github-fe126-cp1-prd.iad.github.net">GitHub</span>, Inc.</li>
        <li><a href="https://github.com/site/terms" data-ga-click="Footer, go to terms, text:terms">Terms</a></li>
        <li><a href="https://github.com/site/privacy" data-ga-click="Footer, go to privacy, text:privacy">Privacy</a></li>
        <li><a href="https://github.com/security" data-ga-click="Footer, go to security, text:security">Security</a></li>
        <li><a href="https://github.com/contact" data-ga-click="Footer, go to contact, text:contact">Contact</a></li>
        <li><a href="https://help.github.com" data-ga-click="Footer, go to help, text:help">Help</a></li>
    </ul>
  </div>
</div>



    
    
    

    <div id="ajax-error-message" class="flash flash-error">
      <span class="octicon octicon-alert"></span>
      <button type="button" class="flash-close js-flash-close js-ajax-error-dismiss" aria-label="Dismiss error">
        <span class="octicon octicon-x"></span>
      </button>
      Something went wrong with that request. Please try again.
    </div>


      <script crossorigin="anonymous" src="https://assets-cdn.github.com/assets/frameworks-161cec7f4cb9a187f66c2ecb8f59c9098e2ad64458af917ec72f9d61d6b8089b.js"></script>
      <script async="async" crossorigin="anonymous" src="https://assets-cdn.github.com/assets/github-24722416f04d8ce310fe98a13efb73e9dd992b51c28b23f229510786915b5f19.js"></script>
      
      
    <div class="js-stale-session-flash stale-session-flash flash flash-warn flash-banner hidden">
      <span class="octicon octicon-alert"></span>
      <span class="signed-in-tab-flash">You signed in with another tab or window. <a href="">Reload</a> to refresh your session.</span>
      <span class="signed-out-tab-flash">You signed out in another tab or window. <a href="">Reload</a> to refresh your session.</span>
    </div>
  </body>
</html>


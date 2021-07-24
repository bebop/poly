remote= curl https://www.uniprot.org/docs/uniprot.xsd | md5sum | cut -d " " -f1 -
local= md5sum uniprot.xsd | cut -d " " -f1 -

if [$remote = $local]
then
	echo All good.
else
	echo OH CRAP IT CHANGED.
fi

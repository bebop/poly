./xsdgen -pkg uniprot uniprot.xsd

sed '/.*Marshal.*/,/^}$/d' xml.go | sed '/.*StatusType) UnmarshalXML.*/,/^}$/d' - | sed '/.*_marshalTime.*/,/^}$/d' - > xml_t.go && mv xml_t.go xml.go

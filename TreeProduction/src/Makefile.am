# One entry point builds the ROOT-only eventTree library and, unless
# configured with --disable-producer, the producer linked to it. Explicit
# source lists; no local-test dependency. Package-qualified headers keep the
# bundled builder separate from CaloReco.
lib_LTLIBRARIES = libEventTree.la
if BUILD_PRODUCER
lib_LTLIBRARIES += libPhotonJetTree.la
endif

AM_CPPFLAGS = -I$(builddir) -I$(srcdir)
# Preserve the staged binary32 arithmetic used by the reference production.
# Fused multiply-add can move a saved shower moment by one float ULP.
AM_CXXFLAGS = -ffp-contract=off

# Saved-product operations: ROOT and the header-only JSON library only.
libEventTree_la_SOURCES = EventTree.cc
nodist_libEventTree_la_SOURCES = EventTreeSchema.h
libEventTree_la_LIBADD = $(ROOT_LIBS)
libEventTree_la_LDFLAGS = -no-undefined -version-info 0:0:0

libPhotonJetTree_la_SOURCES = \
  PhotonJetTree.cc \
  PhotonClusterBuilder.cc \
  internal/Production.cc \
  internal/Setup.cc \
  internal/Event.cc \
  internal/Objects.cc \
  internal/Truth.cc \
  internal/Output.cc
libPhotonJetTree_la_LIBADD = libEventTree.la $(PRODUCER_LIBS)
libPhotonJetTree_la_LDFLAGS = -no-undefined -version-info 0:0:0

productionincludedir = $(includedir)/TreeProduction
productioninclude_HEADERS = EventTree.h
internalincludedir = $(productionincludedir)/internal
internalinclude_HEADERS =
if BUILD_PRODUCER
productioninclude_HEADERS += PhotonJetTree.h PhotonClusterBuilder.h
internalinclude_HEADERS += internal/Production.h internal/Types.h internal/Configuration.h
endif
EXTRA_DIST = autogen.sh ../config/eventTree_v1.yaml

# The one schema file expanded into the branch table compiled into
# libEventTree. Schema::load refuses a schema file with another digest.
schema_file = $(srcdir)/../config/eventTree_v1.yaml
EventTreeSchema.h: $(schema_file)
	$(AM_V_GEN)$(PYTHON) -c 'import hashlib,json,sys,yaml; raw=open(sys.argv[1],"rb").read(); doc=yaml.safe_load(raw); body={k:v for k,v in doc.items() if k!="sources"}; rows=[(n,(n.split("_",1)[0] if "<RR>" in i["group"] else i["group"]),i["type"],i["title"],(i.get("unit","") if s=="branches" else ""),int(s=="reserved")) for s in ("branches","reserved") for i in doc[s] for n in ([i["name"].replace("<RR>",r) for r in ("02","03","04")] if "<RR>" in i["name"] else [i["name"]])]; assert all(x.isascii() and x.isprintable() for r in rows for x in r[:5]), "schema text must be printable ASCII"; print("// Expanded from eventTree_v1.yaml by the build; do not edit."); print("#define EVENTTREE_SCHEMA_SHA256 \"%s\"" % hashlib.sha256(raw).hexdigest()); print("#define EVENTTREE_CONTRACT_SHA256 \"%s\"" % hashlib.sha256(json.dumps(body,sort_keys=True,separators=(",",":")).encode()).hexdigest()); print("const SchemaRow kSchemaRows[] = {"); [print("{%s,%s}," % (",".join(json.dumps(x) for x in r[:5]), "true" if r[5] else "false")) for r in rows]; print("};")' "$(schema_file)" > $@.tmp
	mv $@.tmp $@

# Regenerate when a compiled source, installed header, build recipe or the
# schema changes. The digest describes the actual compilation inputs.
build_identity_sources = $(libEventTree_la_SOURCES) $(libPhotonJetTree_la_SOURCES) EventTree.h PhotonJetTree.h PhotonClusterBuilder.h internal/Production.h internal/Types.h internal/Configuration.h Makefile.am configure.ac autogen.sh ../config/eventTree_v1.yaml
BUILT_SOURCES = EventTreeSchema.h
CLEANFILES = EventTreeSchema.h EventTreeSchema.h.tmp
if BUILD_PRODUCER
BUILT_SOURCES += BuildIdentity.h
CLEANFILES += BuildIdentity.h BuildIdentity.h.tmp
endif
BuildIdentity.h: $(build_identity_sources)
	$(AM_V_GEN)$(PYTHON) -c 'import hashlib,pathlib,sys; r=pathlib.Path(sys.argv[1]); h=hashlib.sha256(); [(h.update(n.encode()+b"\0"),h.update((r/n).read_bytes()),h.update(b"\0")) for n in sorted(sys.argv[2:])]; print("#define PHOTONJET_SOURCE_SHA256 \""+h.hexdigest()+"\"")' "$(srcdir)" $(build_identity_sources) > $@.tmp
	mv $@.tmp $@
